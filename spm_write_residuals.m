function VRes = spm_write_residuals(SPM,Ic)
% Write residual images
% FORMAT Vres = spm_write_residuals(SPM,Ic)
% SPM    - structure containing generic analysis details
% Ic     - contrast index used to adjust data (0:   no adjustment)
%                                             (NaN: adjust for everything) 
%
% VRes   - struct array of residual image handles
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_write_residuals.m 6656 2015-12-24 16:49:52Z guillaume $


%-Get SPM.mat
%--------------------------------------------------------------------------
 if ~nargin || isempty(SPM)
     [SPM,sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
     if ~sts, VRes = ''; return; end
 end

if ~isstruct(SPM)
    swd = spm_file(SPM,'fpath');
    try
        load(fullfile(swd,'SPM.mat'));
        SPM.swd = swd;
    catch
        error(['Cannot read ' fullfile(swd,'SPM.mat')]);
    end
end

try, SPM.swd; catch, SPM.swd = pwd; end
cwd = pwd; cd(SPM.swd);

%-Get contrast used to adjust data
%--------------------------------------------------------------------------
if nargin<2 || isempty(Ic)
    q(1)   = 0;
    Con    = {'<don''t adjust>'};
    q(2)   = NaN;
    Con{2} = '<adjust for everything>';
    for i = 1:length(SPM.xCon)
        if strcmp(SPM.xCon(i).STAT,'F')
            q(end + 1) = i;
            Con{end + 1} = SPM.xCon(i).name;
        end
    end
    i  = spm_input('adjust data for (select contrast)','!+1','m',Con);
    Ic = q(i);
end


%-Compute and write residuals
%==========================================================================
spm('Pointer','Watch')
M   = SPM.xY.VY(1).mat;
DIM = SPM.xY.VY(1).dim(1:min(numel(SPM.xY.VY(1).dim),3));
[nScan, nBeta] = size(SPM.xX.X);

if spm_mesh_detect(SPM.xY.VY)
    file_ext = '.gii';
else
    file_ext = spm_file_ext;
end

%-Initialise residual images
%--------------------------------------------------------------------------
VRes(1:nScan) = deal(struct(...
    'fname',    [],...
    'dim',      DIM,...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'Residuals'));
for i = 1:nScan
    VRes(i).fname   = [sprintf('Res_%04d', i) file_ext];
    VRes(i).descrip = sprintf('Residuals (%04d)', i);
end
VRes = spm_data_hdr_write(VRes);

%-Loop over chunks
%--------------------------------------------------------------------------
chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

spm_progress_bar('Init',nbchunks,'Writing residuals','Chunks');

for i=1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;
    
    %-Get mask
    %----------------------------------------------------------------------
    m = spm_data_read(SPM.VM,chunk) > 0;
    m = m(:)';
    
    %-Get raw data, whiten and filter
    %----------------------------------------------------------------------
    y = zeros(nScan,numel(chunk));
    for j=1:nScan
        y(j,:) = spm_data_read(SPM.xY.VY(j),chunk);
    end
    y(:,~m) = [];
    
    y = spm_filter(SPM.xX.K,SPM.xX.W*y);
    
    if Ic ~= 0
        
        %-Parameter estimates: beta = xX.pKX*xX.K*y
        %------------------------------------------------------------------
        beta = zeros(nBeta,numel(chunk));
        for j=1:nBeta
            beta(j,:) = spm_data_read(SPM.Vbeta(j),chunk);
        end
        beta(:,~m) = [];
        
        %-Subtract Y0 = XO*beta,  Y = Yc + Y0 + e
        %------------------------------------------------------------------
        if ~isnan(Ic)
            y = y - spm_FcUtil('Y0',SPM.xCon(Ic),SPM.xX.xKXs,beta);
        else
            y = y - SPM.xX.xKXs.X * beta;
        end
        
    end
    
    %-Write residuals
    %----------------------------------------------------------------------
    yy = NaN(numel(chunk),1);
    for j=1:nScan
        yy(m)   = y(j,:);
        VRes(j) = spm_data_write(VRes(j), yy, chunk); 
    end
    
    spm_progress_bar('Set',i)
end

cd(cwd);
spm_progress_bar('Clear')
spm('Pointer','Arrow')

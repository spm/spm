function VRes = spm_write_residuals(SPM,Ic)
% Write residual images
% FORMAT Vres = spm_write_residuals(SPM,Ic)
% SPM    - structure containing generic analysis details
% Ic     - contrast index used to adjust data (0:   no adjustment)
%                                             (NaN: adjust for everything) 
%
% Vres   - struct array of residual image handles
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_write_residuals.m 5549 2013-06-12 12:41:18Z gareth $


%-Get SPM.mat
 %--------------------------------------------------------------------------
 if ~nargin || isempty(SPM)
     [SPM,sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
     if ~sts, Vres = ''; return; end
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
N   = numel(SPM.xY.VY);
M   = SPM.xY.VY(1).mat;
DIM = SPM.xY.VY(1).dim(1:3)';

%-Initialise residual images
%--------------------------------------------------------------------------
VRes(1:N) = deal(struct(...
    'fname',    [],...
    'dim',      DIM',...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'Residuals'));
for i = 1:N
    VRes(i).fname   = [sprintf('Res_%04d', i) spm_file_ext];
    VRes(i).descrip = sprintf('Residuals (%04d)', i);
end
VRes = spm_create_vol(VRes);

%-Loop over slices
%--------------------------------------------------------------------------
spm_progress_bar('Init',DIM(3),'Slices');
for i=1:DIM(3)
    
    %-Get mask
    %----------------------------------------------------------------------
    m = spm_slice_vol(SPM.VM,spm_matrix([0 0 i]),DIM(1:2),0) > 0;
    m = m(:)';
    
    %-Get raw data, whiten and filter
    %----------------------------------------------------------------------
    y = zeros(N,prod(DIM(1:2)));
    for j=1:N
        tmp = spm_slice_vol(SPM.xY.VY(j),spm_matrix([0 0 i]),DIM(1:2),0);
        y(j,:) = tmp(:)';
    end
    y(:,~m) = [];
    
    y = spm_filter(SPM.xX.K,SPM.xX.W*y);
    
    if Ic ~= 0
        
        %-Parameter estimates: beta = xX.pKX*xX.K*y
        %------------------------------------------------------------------
        beta = zeros(numel(SPM.Vbeta),prod(DIM(1:2)));
        for j=1:numel(SPM.Vbeta)
            tmp = spm_slice_vol(SPM.Vbeta(j),spm_matrix([0 0 i]),DIM(1:2),0);
            beta(j,:) = tmp(:)';
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
    yy = NaN(DIM(1:2)');
    for j=1:N
        yy(m) = y(j,:);
        VRes(j) = spm_write_plane(VRes(j),yy,i);
    end
    
    spm_progress_bar('Set',i)
end

cd(cwd);
spm_progress_bar('Clear')
spm('Pointer','Arrow')

function [xCon,SPM]= spm_vb_x2(SPM,XYZ,xCon,ic)
% Compute and write Chi^2 image
% FORMAT [xCon,SPM]= spm_vb_x2(SPM,XYZ,xCon,ic)
%
% SPM  - SPM data structure
% XYZ  - voxel list
% xCon - contrast info
% ic   - contrast number
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_x2.m 2451 2008-11-10 16:20:32Z lee $

% Get approximate posterior covariance for ic
% using Taylor-series approximation
        
%-Get number of sessions
%-----------------------------------------------------------------------
nsess = length(SPM.nscan); %length(SPM.Sess);

%- Compound Contrast
%-----------------------------------------------------------------------
c  = xCon(ic).c;
kc = size(c,2);

%-Get posterior beta's
%-----------------------------------------------------------------------
Nk = size(SPM.xX.X,2);

for k=1:Nk,
    beta(k,:) = spm_get_data(SPM.VCbeta(k),XYZ);
end

%-Get posterior SD beta's
%-----------------------------------------------------------------------
Nk = size(SPM.xX.X,2);

for k=1:Nk,
    sd_beta(k,:) = spm_get_data(SPM.VPsd(k),XYZ);
end

%-Get AR coefficients
%-----------------------------------------------------------------------
for s=1:nsess,
    for p=1:SPM.PPM.AR_P,
        Sess(s).a(p,:) = spm_get_data(SPM.PPM.Sess(s).VAR(p),XYZ);
    end
end

%-Get noise SD
%-----------------------------------------------------------------------
for s=1:nsess,
    Sess(s).lambda = spm_get_data(SPM.PPM.Sess(s).VHp,XYZ);
end

%-Loop over voxels
%=======================================================================
Nvoxels = size(XYZ,2);
D       = repmat(NaN,reshape(SPM.xVol.DIM(1:3),1,[]));

spm_progress_bar('Init',100,'Estimating posterior contrast variance','');

for v=1:Nvoxels,
    %-Which block are we in ?
    %-------------------------------------------------------------------
    block_index = SPM.xVol.labels(1,v);
    
    V = zeros(kc,kc);
    m = zeros(kc,1);
    d = 0;
    for s=1:nsess,
        
        %-Reconstruct approximation to voxel wise correlation matrix
        %---------------------------------------------------------------
        R = SPM.PPM.Sess(s).block(block_index).mean.R;
        if SPM.PPM.AR_P > 0
            dh = Sess(s).a(:,v)'-SPM.PPM.Sess(s).block(block_index).mean.a;
            dh = [dh Sess(s).lambda(v)-SPM.PPM.Sess(s).block(block_index).mean.lambda];
            for i=1:length(dh),
                R = R + SPM.PPM.Sess(s).block(block_index).mean.dR(:,:,i) * dh(i);
            end 
        end
        %-Get indexes of regressors specific to this session
        %---------------------------------------------------------------
        scol           = SPM.Sess(s).col; 
        mean_col_index = SPM.Sess(nsess).col(end) + s;
        scol           = [scol mean_col_index];
        
        %-Reconstruct approximation to voxel wise covariance matrix
        %---------------------------------------------------------------
        Sigma_post = (sd_beta(scol,v) * sd_beta(scol,v)') .* R;
        
        % Get component of contrast covariance specific to this session
        %---------------------------------------------------------------
        CC = c(scol,:);
        V  = V + CC' * Sigma_post * CC; 
        
        % Get posterior mean contrast vector
        m = m + CC'*beta(scol,v);
        
        % Get Chi^2 value
        d = d + m'*inv(V)*m;
    end
    
    D(XYZ(1,v),XYZ(2,v),XYZ(3,v)) = d;
    if rem(v,100)==0
        % update progress bar every 100th voxel
        spm_progress_bar('Set',100*v/Nvoxels);
    end
    
end

xCon(ic).eidf=rank(V);

spm_progress_bar('Clear');   

% Create handle
%-----------------------------------------------------------------------
Vhandle = struct(...
    'fname',  sprintf('x2_%04d.img',ic),...
    'dim',    SPM.xVol.DIM',...
    'dt',     [spm_type('float32') spm_platform('bigend')],... 
    'mat',    SPM.xVol.M,...
    'pinfo',  [1,0,0]',...
    'descrip',sprintf('Chi^2 stat for Bayes multivar con %d: %s',ic,xCon(ic).name));

%-Write image
%-----------------------------------------------------------------------
Vhandle = spm_create_vol(Vhandle);
Vhandle = spm_write_vol(Vhandle,D);

xCon(ic).Vcon = Vhandle;

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),...
    sprintf('...written %s',spm_str_manip(Vhandle.fname,'t')));            %-#

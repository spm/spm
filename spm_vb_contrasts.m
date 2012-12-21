function SPM = spm_vb_contrasts(SPM,XYZ,xCon,ic)
% Compute and write posterior standard deviation image of given contrast
% FORMAT SPM = spm_vb_contrasts(SPM,XYZ,xCon,ic)
%
% SPM   - SPM data structure
% XYZ   - voxel list
% xCon  - contrast info
% ic    - contrast number
%
% Get approximate posterior covariance for given contrast ic using 
% Taylor-series approximation
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_vb_contrasts.m 5160 2012-12-21 16:58:38Z guillaume $

        
%-Get number of sessions, beta's and voxels
%--------------------------------------------------------------------------
nsess   = length(SPM.Sess);
nbeta   = size(SPM.xX.X,2);
nvox    = size(XYZ,2);

%-Get contrast vector
%--------------------------------------------------------------------------
c       = xCon(ic).c;

%-Get posterior SD beta's
%--------------------------------------------------------------------------
sd_beta = zeros(nbeta,nvox);
for k=1:nbeta
    sd_beta(k,:) = spm_get_data(SPM.VPsd(k),XYZ);
end

%-Get AR coefficients
%--------------------------------------------------------------------------
[Sess(1:nsess)] = deal(struct('a',zeros(SPM.PPM.AR_P,nvox),'lambda',[]));
for s=1:nsess
    for p=1:SPM.PPM.AR_P
        Sess(s).a(p,:) = spm_get_data(SPM.PPM.Sess(s).VAR(p),XYZ);
    end
end

%-Get noise SD
%--------------------------------------------------------------------------
for s=1:nsess
    Sess(s).lambda = spm_get_data(SPM.PPM.Sess(s).VHp,XYZ);
end

%-Loop over voxels
%==========================================================================
Y       = NaN(reshape(SPM.xVol.DIM(1:3),1,[]));

spm_progress_bar('Init',100,'Estimating posterior contrast variance','');

for v=1:nvox
    %-Which block are we in ?
    %----------------------------------------------------------------------
    block_index = SPM.xVol.labels(1,v);
    
    y = 0;
    for s=1:nsess
        
        %-Reconstruct approximation to voxel wise correlation matrix
        %------------------------------------------------------------------
        R = SPM.PPM.Sess(s).block(block_index).mean.R;
        if SPM.PPM.AR_P > 0
            dh = Sess(s).a(:,v)'-SPM.PPM.Sess(s).block(block_index).mean.a;
            dh = [dh Sess(s).lambda(v)-SPM.PPM.Sess(s).block(block_index).mean.lambda];
            for i=1:length(dh),
                R = R + SPM.PPM.Sess(s).block(block_index).mean.dR(:,:,i) * dh(i);
            end 
        end
        
        %-Get indexes of regressors specific to this session
        %------------------------------------------------------------------
        scol           = SPM.Sess(s).col; 
        mean_col_index = SPM.Sess(nsess).col(end) + s;
        scol           = [scol mean_col_index];
        
        % Revert from voxel-specific to slice-specific R if Taylor-series
        % approximation produces invalid R
        if max(max(abs(R))) > 1
            R=SPM.PPM.Sess(s).block(block_index).mean.R;
        end
        
        %-Reconstruct approximation to voxel wise covariance matrix
        %------------------------------------------------------------------
        Sigma_post = (sd_beta(scol,v) * sd_beta(scol,v)') .* R;
        
        % Get component of contrast variance specific to this session
        %------------------------------------------------------------------
        CC = c(scol);
        y  = y + CC' * Sigma_post * CC; 
        
    end
    
    Y(XYZ(1,v),XYZ(2,v),XYZ(3,v)) = sqrt(y);
    if rem(v,100)==0
        % update progress bar every 100th voxel
        spm_progress_bar('Set',100*v/nvox);
    end
    
end

spm_progress_bar('Clear');   

% Create handle
%--------------------------------------------------------------------------
V = struct(...
    'fname',  [sprintf('con_sd_%04d',ic) spm_file_ext],...
    'dim',    SPM.xVol.DIM',...
    'dt',     [spm_type('float32') spm_platform('bigend')],... 
    'mat',    SPM.xVol.M,...
    'pinfo',  [1,0,0]',...
    'descrip',sprintf('PPM contrast SD - %d: %s',ic,xCon(ic).name));

%-Write image
%--------------------------------------------------------------------------
V = spm_create_vol(V);
V = spm_write_vol(V,Y);

SPM.PPM.Vcon_sd(ic) = V;

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),...
    sprintf('...written %s',spm_file(V.fname,'filename')));             %-#

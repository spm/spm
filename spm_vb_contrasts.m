function [SPM]= spm_vb_contrasts(SPM,XYZ,xCon,ic)
% Compute and write posterior SD image
% FORMAT [SPM]= spm_vb_contrasts(SPM,XYZ,xCon,ic)
%
% SPM   SPM data structure
% XYZ   voxel list
% xCon  Contrast info
% ic    contrast number
%
% xCon  Updated contrast info
%
% Will Penny $Id$


% Get approximate posterior covariance for ic
% using Taylor-series approximation
        
% Get number of sessions
nsess=length(SPM.Sess);

% Contrast
c     = xCon(ic).c;

% Get posterior SD beta's
Nk=size(SPM.xX.X,2);

for k=1:Nk,
    sd_beta(k,:) = spm_get_data(SPM.VPsd(k),XYZ);
end


% Get AR coefficients
for s=1:nsess,
    for p=1:SPM.PPM.AR_P
        Sess(s).a(p,:) = spm_get_data(SPM.PPM.Sess(s).VAR(p),XYZ);
    end
end

% Get noise SD
for s=1:nsess,
    Sess(s).lambda = spm_get_data(SPM.PPM.Sess(s).VHp,XYZ);
end

% Loop over voxels
Nvoxels=size(XYZ,2);
xdim=SPM.xVol.DIM(1);
ydim=SPM.xVol.DIM(2);
zdim=SPM.xVol.DIM(3);
Y  = NaN*ones(xdim,ydim,zdim);
spm_progress_bar('Init',100,'Estimating posterior contrast variance','');
for v=1:Nvoxels,
    % Which slice are we in ?
    slice_index=XYZ(3,v);
    
    y=0;
    for s=1:nsess,
        
        % Reconstruct approximation to voxel wise correlation matrix
        R=SPM.PPM.Sess(s).slice(slice_index).mean.R;
        dh=Sess(s).a(:,v)'-SPM.PPM.Sess(s).slice(slice_index).mean.a;
        dh=[dh Sess(s).lambda(v)-SPM.PPM.Sess(s).slice(slice_index).mean.lambda];
        for i=1:length(dh),
            R=R+SPM.PPM.Sess(s).slice(slice_index).mean.dR(:,:,i)*dh(i);
        end 
        
        % Get indexes of regressors specific to this session
        scol=SPM.Sess(s).col; 
        mean_col_index=SPM.Sess(nsess).col(end)+s;
        scol=[scol mean_col_index];
        
        % Reconstruct approximation to voxel wise covariance matrix
        Sigma_post = (sd_beta(scol,v)*sd_beta(scol,v)').*R;
        
        % Get component of contrast variance specific to this session
        CC=c(scol);
        y=y+CC'*Sigma_post*CC; 
    end
    Y(XYZ(1,v),XYZ(2,v),XYZ(3,v)) = sqrt(y);
    spm_progress_bar('Set',100*v/Nvoxels);
end
spm_progress_bar('Clear');   

% Create handle
V = struct(...
    'fname',  sprintf('con_sd_%04d.img',ic),...
    'dim',    SPM.xVol.DIM',...
    'dt',     [16, spm_platform('bigend')],...
    'mat',    SPM.xVol.M,...
    'pinfo',  [1,0,0]',...
    'descrip',sprintf('PPM contrast SD - %d: %s',ic,xCon(ic).name));

%-Write image
%-----------------------------------------------------------
V = spm_create_vol(V,'noopen');
V  = spm_write_vol(V,Y);
V = spm_close_vol(V);

SPM.PPM.Vcon_sd(ic)=V;

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
    '...written %s',spm_str_manip(V.fname,'t')))%-#

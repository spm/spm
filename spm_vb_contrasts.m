function [SPM]= spm_vb_contrasts(SPM,XYZ,xCon,ic)
% Compute and write posterior variance image
% FORMAT [SPM]= spm_vb_contrasts(SPM,XYZ,xCon,ic)
%
% SPM   SPM data structure
% XYZ   voxel list
% xCon  Contrast info
% ic    contrast number
%
% xCon  Updated contrast info

% Get approximate posterior covariance for ic
% using Taylor-series approximation
        
% Contrast
c     = xCon(ic).c;

% Get posterior SD beta's
Nk=size(SPM.xX.X,2);
for k=1:Nk,
    sd_beta(k,:) = spm_get_data(SPM.VPsd(k),XYZ);
end

% Get AR coefficients
for p=1:SPM.PPM.AR_P
    a(p,:) = spm_get_data(SPM.VAR(p),XYZ);
end

% Get noise SD
lambda = spm_get_data(SPM.VHp,XYZ);

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
    
    % Reconstruct approximation to voxel wise correlation matrix
    R=SPM.PPM.slice(slice_index).mean.R;
    dh=a(:,v)'-SPM.PPM.slice(slice_index).mean.a;
    dh=[dh lambda(v)-SPM.PPM.slice(slice_index).mean.lambda];
    for i=1:length(dh),
        R=R+SPM.PPM.slice(slice_index).mean.dR(:,:,i)*dh(i);
    end 
    % Reconstruct approximation to voxel wise covariance matrix
    Sigma_post = (sd_beta(:,v)*sd_beta(:,v)').*R;
    Y(XYZ(1,v),XYZ(2,v),XYZ(3,v)) = c'*Sigma_post*c;
    spm_progress_bar('Set',100*v/Nvoxels);
end
spm_progress_bar('Clear');   

% Create handle
V = struct(...
    'fname',  sprintf('con_var_%04d.img',ic),...
    'dim',    [SPM.xVol.DIM',16],...
    'mat',    SPM.xVol.M,...
    'pinfo',  [1,0,0]',...
    'descrip',sprintf('PPM contrast variance - %d: %s',ic,xCon(ic).name));

%-Write image
%-----------------------------------------------------------
V = spm_create_vol(V,'noopen');
V  = spm_write_vol(V,Y);
V = spm_close_vol(V);

%xCon(ic).Vcon_var=V;
SPM.PPM.Vcon_var(ic)=V;

fprintf('%s%30s\n',sprintf('\b')*ones(1,30),sprintf(...
    '...written %s',spm_str_manip(V.fname,'t')))%-#
function [noise,mu_val] = spm_noise_estimate(Scans)
% Estimate avarage noise from a series of images
% FORMAT noise = spm_noise_estimate(Scans)
% Scans  - nifti structures or filenames of images
% noise  - standard deviation estimate
% mu_val - expectation of more intense Rician
% _______________________________________________________________________
%  Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% $Id: spm_noise_estimate.m 7460 2018-10-29 15:55:12Z john $

if ~isa(Scans,'nifti'), Scans = nifti(Scans); end

noise  = zeros(numel(Scans),1);
mu_val = zeros(numel(Scans),1);
for i=1:numel(Scans)
    Nii = Scans(i);
    f   = Nii.dat(:,:,:);
    if spm_type(Nii.dat.dtype(1:(end-3)),'intt')
        f(f==max(f(:))) = 0;
        x      = 0:Nii.dat.scl_slope:max(f(:));
        [h,x]  = hist(f(f~=0),x);
    else
        x      = (0:1023)*(max(f(:))/1023);
        f(f==max(f(:))) = 0;
        [h,x]  = hist(f(f~=0 & isfinite(f)),x);
    end
    [mg,nu,sd] = spm_rice_mixture(double(h(:)),double(x(:)),2);
    noise(i)   = min(sd);

    if nargout>=2
        x          = -nu.^2./(2*sd.^2);
        msk        = x>-20;
        Laguerre   = exp(x(msk)/2).*((1-x(msk)).*besseli(0,-x(msk)/2)-x(msk).*besseli(1,-x(msk)/2));
        Ey         = zeros(size(sd));
        Ey( msk)   = sqrt(pi*sd(msk).^2/2).*Laguerre;
        Ey(~msk)   = nu(~msk);
        mu_val(i)  = max(Ey);
    end
end


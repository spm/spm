function noise = spm_noise_estimate(Scans)
% Estimate avarage noise from a series of images
% FORMAT noise = spm_noise_estimate(Scans)
% Scans - nifti structures or filenames of images
% noise - standard deviation estimate
% _______________________________________________________________________
%  Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% $Id: spm_noise_estimate.m 4776 2012-07-02 20:33:35Z john $

if ~isa(Scans,'nifti'), Scans = nifti(Scans); end

noise = zeros(numel(Scans),1);
for i=1:numel(Scans),
    Nii = Scans(i);
    f   = Nii.dat(:,:,:);
    if spm_type(Nii.dat.dtype(1:(end-3)),'intt'),
        f(f==max(f(:))) = 0;
        x      = 0:Nii.dat.scl_slope:max(f(:));
        [h,x]  = hist(f(f~=0),x);
    else
        x      = (0:1023)*(max(f(:))/1023);
        f(f==max(f(:))) = 0;
        [h,x]  = hist(f(f~=0 & isfinite(f)),x);
    end
    [mg,nu,sd] = spm_rice_mixture(h(:),x(:),2);
    noise(i)   = min(sd);
end


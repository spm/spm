function [noise,mu_val,info] = spm_noise_estimate(Scans,K)
% Estimate average noise from a series of images
% FORMAT [noise,mu_val] = spm_noise_estimate(Scans)
% Scans  - nifti structures or filenames of images
% K      - Number of Rician mixture components
% noise  - standard deviation estimate
% mu_val - expectation of more intense Rician
% info   - This struct can be used for plotting the fit as:
%              plot(info.x(:),info.p,'--',info.x(:), ...
%                   info.h/sum(info.h)/info.md,'b.', ...
%                   info.x(:),info.sp,'r');
% _________________________________________________________________________
%  Copyright (C) 2012-2019 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_noise_estimate.m 7599 2019-05-30 13:50:41Z mikael $

if ~isa(Scans,'nifti'), Scans = nifti(Scans); end

if nargin < 2, K = 2; end

noise  = zeros(numel(Scans),1);
mu_val = zeros(numel(Scans),1);
info   = struct('x',[],'h',[],'p',[],'sp',[],'md',[]);
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
    [mg,nu,sd,info(i)] = spm_rice_mixture(double(h(:)),double(x(:)),K);
    noise(i)           = min(sd);

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

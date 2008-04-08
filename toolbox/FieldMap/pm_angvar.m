function angvar = pm_angvar(cmap)
% Estimates the (voxelwise) variance of the angle
% estimated from the complex map cmap.
% FORMAT: angvar = pm_angvar(cmap)
%
% Input:
% cmap     : Complex-valued MR intensity image. When used to
%            estimate the variance of a delta_phi map estimated
%            from two measurements with different echo-time this
%            should be the image with the longer echo-time.
%
% Output:
% angvar   : Map with an estimate of the variance of a phasemap
%            estimated using cmap as one of its constituents.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_angvar.m 1317 2008-04-08 16:16:38Z chloe $

% Get mutual histogram of Re and Im part of all voxels.
%
maxr = max(real(cmap(:)));
minr = min(real(cmap(:)));
if maxr > abs(minr) minr = -maxr; else maxr = -minr; end
maxi = max(imag(cmap(:)));
mini = min(imag(cmap(:)));
if maxi > abs(mini) mini = -maxi; else maxi = -mini; end
hs = 255;


MH = sparse(round((hs-1)*((imag(cmap(:))-mini)/(maxi-mini))+1),...
            round((hs-1)*((real(cmap(:))-minr)/(maxr-minr))+1),...
            ones(prod(size(cmap)),1),hs,hs);

%
% Estimate variance of noise-peak, utilising the
% fact that the tissue peak has been spread out over 2pi,
% and hence contribute very little.
%

pdf = MH(ceil(hs/2),:)/sum(MH(ceil(hs/2),:));
variance(1) = full(sum(pdf.*((([1:hs]-ceil(hs/2))*(maxr-minr)/hs).^2)));
pdf = MH(:,ceil(hs/2))/sum(MH(:,ceil(hs/2)));
variance(2) = full(sum(pdf.*((([1:hs]'-ceil(hs/2))*(maxi-mini)/hs).^2)));
stdev = sqrt(mean(variance));

%
% Use the estimate of variance to do a crude simulation
% of the uncertainty of the angle measurement.
%

randn('state',0);       % Make sure results are reproducible.
avals = linspace(0,max(abs(cmap(:))),100);
sim = repmat(avals/sqrt(2) + sqrt(-1)*avals/sqrt(2),1000,1);
sim1 = sim + complex(stdev*randn(1000,100),stdev*randn(1000,100));
sim2 = sim + complex(stdev*randn(1000,100),stdev*randn(1000,100));
sim = angle(exp(sqrt(-1)*angle(sim1))./exp(sqrt(-1)*angle(sim2)));
sim(find(sim(:) < -pi)) = sim(find(sim(:) < -pi)) + 2*pi;
sim(find(sim(:) > pi)) = sim(find(sim(:) > pi)) - 2*pi;
variance = var(sim);
angvar = reshape(interp1(avals,variance,abs(cmap(:)),'linear'),size(cmap));
 
%
% Smooth ever so little
%

%spm_smooth(angvar,angvar,1.5);
%angvar(1,:,:) = 1.2*angvar(1,:,:); angvar(end,:,:) = 1.2*angvar(end,:,:);
%angvar(:,1,:) = 1.2*angvar(:,1,:); angvar(:,end,:) = 1.2*angvar(:,end,:);

return;





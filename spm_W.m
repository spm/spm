function [W,FWHM] = spm_W(V)
% returns the smoothness estimator for a 3-D memory mapped image
% FORMAT [W FWHM] = spm_W(V);
% V    - memory mapped volume
%
% W    - {1 x 3} vector of smoothness estimates in x y and z {voxels}
% FWHM - equivalent estimator expressed as FWHM {mm}
%____________________________________________________________________________
%
% spm_W returns the smoothness estimator for a 3-D memory mapped image
% in voxels using the variance of the process and its first partial
% derivatives.
%
%__________________________________________________________________________
% %W% %E%

%--------------------------------------------------------------------------
D     = zeros(V(1,1)*V(2,1),1);			% dummy matrix for smoothness
sx    = zeros(1,2);				% smoothness estimators {x}
sy    = zeros(1,2);				% smoothness estimators {y}
sz    = zeros(1,2);				% smoothness estimators {z}

for i = 1:V(3,1)

	% sums of squares of SPM{t} and spatial derivatives
        %-------------------------------------------------------------------
	d       = spm_slice_vol(V,spm_matrix([0 0 i]),[V(1,1) V(2,1)],0);
	d       = d - mean(d(:));
	dz      = D - d(:);
	[dy dx] = gradient(d);
	Y       = ~d;
	Y       = ~(Y | abs(gradient(Y))); 
	Y       = Y(:);
	sx(1)   = sx(1) + sum( d(Y).^2);
	sx(2)   = sx(2) + sum(dx(Y).^2);
	sy(1)   = sy(1) + sum( d(Y).^2);
	sy(2)   = sy(2) + sum(dy(Y).^2);
	sz(1)   = sz(1) + sum( d(D & d(:)).^2);
	sz(2)   = sz(2) + sum(dz(D & d(:)).^2);
	D       = d(:);
end
W	  = sqrt([sx(:,1)./sx(:,2) sy(:,1)./sy(:,2) sz(:,1)./sz(:,2)]/2);

%---------------------------------------------------------------------------
if V(3,1) == 1;   W = W(1:2);  end			% 2 dimnesional data

FWHM      = sqrt(8*log(2))*W.*V(([1:length(W)] + 3),1)';

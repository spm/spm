function [Transform,Dims,remainder] = spm_snbasis(VG,VF,Affine,param)

% 3D Basis Function Normalization
% FORMAT [Transform,Dims,remainder] = spm_snbasis(VG,VF,Affine,param);
% VG        - Vector of memory mapped templates.
% VF        - Memory mapped image to normalize.
% Affine    - A 4x4 transformation (in voxel space).
% param(1:3)- Number of basis functions in each dimension.
% param(4)  - Number of iterations.
% param(5)  - smoothness of the image.
% Transform - Three column vectors containing the discrete cosine
%             transform of the warps in X, Y & Z.
% Dims      - The dimensions of the transforms.
% remainder - The scaling for each of the templates required
%             to get the best match to the normalized image.
%
% spm_snbasis performs a spatial normalization based upon a 3D
% discrete cosine transform.
% The images should all have an 8mm FWHM smoothness.
%

% %W% John Ashburner FIL %E%

global Alpha Beta T IC0 Var

fwhm = param(5);
sd1  = param(6);

% Number of basis functions for x, y & z
%-----------------------------------------------------------------------
k=param(1:3);
k = max([ k ; 1 1 1]);
k = min([ k ; VG(1) VG(2) VG(3)]);

% Scaling is to improve stability.
%-----------------------------------------------------------------------
stabilise = 64;
basX = spm_dctmtx(VG(1),k(1))*stabilise;
basY = spm_dctmtx(VG(2),k(2))*stabilise;
basZ = spm_dctmtx(VG(3),k(3))*stabilise;

dbasX = spm_dctmtx(VG(1),k(1),'diff')*stabilise;
dbasY = spm_dctmtx(VG(2),k(2),'diff')*stabilise;
dbasZ = spm_dctmtx(VG(3),k(3),'diff')*stabilise;

% Estimate a suitable sparse diagonal inverse covariance matrix for
% the parameters (IC0). These are based on gradients of deformation
% functions, such that the average gradient of deformation fields
% generated from random normal distributions of parameters with
% covariance matrix inv(IC0) is approximately equal to sd1.
%-----------------------------------------------------------------------
kk=(2:k(1))'; dd1 = zeros(k(1),1); dd1(kk) = (pi*(kk-1)/VG(1)).^2;
kk=(2:k(2))'; dd2 = zeros(k(2),1); dd2(kk) = (pi*(kk-1)/VG(2)).^2;
kk=(2:k(3))'; dd3 = zeros(k(3),1); dd3(kk) = (pi*(kk-1)/VG(3)).^2;

IC0 = 	kron(kron(ones(k(3),1),ones(k(2),1)),dd1         ) + ...
	kron(kron(ones(k(3),1),dd2         ),ones(k(1),1)) + ...
	kron(kron(dd3         ,ones(k(2),1)),ones(k(1),1));
IC0 = IC0/(3*(prod(VG(1:3))/(prod(k)-1)));
IC0 = IC0*(sd1.^(-2))*stabilise^6;
IC0 = [IC0 ; IC0 ; IC0 ; zeros(size(VG,2)*4,1)];
IC0 = sparse(1:length(IC0),1:length(IC0),IC0,length(IC0),length(IC0));


% Generate starting estimates.
%-----------------------------------------------------------------------
s1 = 3*prod(k);
s2 = s1 + size(VG,2)*4;
T = zeros(s2,1);
T(s1+(1:4:size(VG,2)*4)) = ones(size(VG,2),1);

for iter=1:param(4)
	fprintf('iteration # %d: ', iter);
	[Alpha,Beta,Var] = spm_brainwarp(VG,VF,Affine,basX,basY,basZ,dbasX,dbasY,dbasZ,T,fwhm);
	if (iter > 0)
		% Parameter estimates biased towards affine.
		%
		% Assume that the problem is linear.
		% The unbiased solution would be:
		% X = T + Alpha\Beta;
		%
		% If we assume that each observation F has the same amount
		% of Gaussian noise, we can approximate the variance of this noise
		% by:
		% Var = \sum((F - A*T)^2)/(n-1) - where n is something like
		% the number of 'resels' in the image F, and A is the design matrix.
		% The formal covariance matrix of the solution X, is:
		% inv(Alpha)*Var
		%
		% The covariance matrix describing the distribution of warps among
		% the population is estimated as inv(IC0).
		% These warps are deviations away from the affine, ie. the origin
		% of the parameter space.
		%
		% The solution we require is the weighted mean of 0 and T + Alpha\Beta,
		% where the weights are the inverses of the covariance matrixes:
		%
		% ie.
		% (inv(C1) + inv(C2))\(inv(C1)*X1 + inv(C2)*X2)
		% where:	X1 = T + Alpha\Beta
		% 		C1 = inv(Alpha)*Var
		% 		X2 = 0
		% 		C2 = inv(IC0)
		%
		% Which simplifies to: (Alpha + IC0*Var)\(Alpha*T + Beta)
		%-----------------------------------------------------------------------
		T = (Alpha + IC0*Var)\(Alpha*T + Beta);
	else
		% The simple straightforward solution.
		T = T + Alpha\Beta;
	end
	disp(Var);
end

% Dimensions and values of the 3D-DCT
%-----------------------------------------------------------------------
Transform = reshape(T(1:s1),prod(k),3)*stabilise.^3;
Dims = [VG(1:3,1)'; k];

% Scaling for each template image.
%-----------------------------------------------------------------------
remainder = T((1:4:size(VG,2)*4) + s1);

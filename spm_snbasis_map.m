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


fwhm = param(5);
sd1  = param(6);

% Number of basis functions for x, y & z
%-----------------------------------------------------------------------
k=param(1:3);
k = max([ k ; 1 1 1]);
k = min([ k ; VG(1) VG(2) VG(3)]);

% Scaling is to improve stability.
%-----------------------------------------------------------------------
stabilise = 200;
basX = spm_dctmtx(VG(1),k(1))*stabilise;
basY = spm_dctmtx(VG(2),k(2))*stabilise;
basZ = spm_dctmtx(VG(3),k(3))*stabilise;

% Guess the Inverse Covariance matrix for brain deformations.
% These are based on gradients of deformation functions
% in the direction of the warp.
% eg diff(basX*(VX1.*randn(size(VX1)))) should have a mean RMS amplitude of sd1;
%-----------------------------------------------------------------------
VX1 = (0.001+(0:k(1)-1))'.^(-1) * sd1 * VG(1)^(3/2) / pi / sqrt(0.001+(k(1)-1))/stabilise;
VY1 = (0.001+(0:k(2)-1))'.^(-1) * sd1 * VG(2)^(3/2) / pi / sqrt(0.001+(k(2)-1))/stabilise;
VZ1 = (0.001+(0:k(3)-1))'.^(-1) * sd1 * VG(3)^(3/2) / pi / sqrt(0.001+(k(3)-1))/stabilise;

% eg basX*(VX2.*randn(size(VX2))) should have a mean RMS amplitude of 1;
% some type of Hanning like wrapper may be used in the future
%-----------------------------------------------------------------------
VX2 = ones(k(1),1) * sqrt(VG(1)/k(1)) / stabilise;
VY2 = ones(k(2),1) * sqrt(VG(2)/k(2)) / stabilise;
VZ2 = ones(k(3),1) * sqrt(VG(3)/k(3)) / stabilise;

% multiply together and store on diagonal of a sparse matrix.
%-----------------------------------------------------------------------
IC0 = [kron(VZ2,kron(VY2,VX1)); kron(VZ2,kron(VY1,VX2)); kron(VZ1,kron(VY2,VX2))].^(-2);
IC0 = sparse(1:length(IC0),1:length(IC0),IC0,length(IC0)+size(VG,2)*4,length(IC0)+size(VG,2)*4);


% Generate starting estimates.
%-----------------------------------------------------------------------
s1 = 3*prod(k);
s2 = s1 + size(VG,2)*4;
T = zeros(s2,1);
T(s1+(1:4:size(VG,2)*4)) = ones(size(VG,2),1);

for iter=1:param(4)
	fprintf('iteration # %d: ', iter);
	[Alpha,Beta,Var] = spm_brainwarp2(VG,VF,Affine,basX,basY,basZ,T,fwhm);
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
	%disp(Var);
end

% Dimensions and values of the 3D-DCT
%-----------------------------------------------------------------------
Transform = reshape(T(1:s1),prod(k),3)*stabilise.^3;
Dims = [VG(1:3,1)'; k];

% Scaling for each template image.
%-----------------------------------------------------------------------
remainder = T((1:4:size(VG,2)*4) + s1);

% 3D Basis Function Normalization
% FORMAT [Transform,Dims,remainder] = spm_snbasis(VG,VF,Affine,param,fwhm);
% VG        - Vector of memory mapped templates.
% VF        - Memory mapped image to normalize.
% Affine    - A 4x4 transformation (in voxel space).
% param(1:3)- Number of basis functions in each dimension.
% param(4)  - Number of iterations.
% fwhm      - smoothness of the image.
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
% This work is still in progress. Future work is likely to be involved
% in constraining the "energy" in the deformation field.

% %W% John Ashburner FIL %E%

function [Transform,Dims,remainder] = spm_snbasis(VG,VF,Affine,param,fwhm)

fwhm = param(5);
sd1  = param(6);

% Number of basis functions for x, y & z
%k1 = floor(VG(1:3,:).*VG(4:6,:)/(2*param(1))+1);
k=param(1:3);
k = max([ k ; 1 1 1]);
k = min([ k ; VG(1) VG(2) VG(3)]);

stabilise = 200;
% Scaling is to try and improve stability.
basX = spm_dctmtx(VG(1),k(1))*stabilise;
basY = spm_dctmtx(VG(2),k(2))*stabilise;
basZ = spm_dctmtx(VG(3),k(3))*stabilise;

% Guess the Inverse Covariance matrix for brain deformations.
% These are based on gradients of deformation functions
% in the direction of the warp.
% eg diff(basX*(VX1.*randn(size(VX1)))) should have a mean RMS amplitude of sd1;
VX1 = (0.001+(0:k(1)-1))'.^(-1) * sd1 * VG(1)^(3/2) / pi / sqrt((k(1)-1))/stabilise;
VY1 = (0.001+(0:k(2)-1))'.^(-1) * sd1 * VG(2)^(3/2) / pi / sqrt((k(2)-1))/stabilise;
VZ1 = (0.001+(0:k(3)-1))'.^(-1) * sd1 * VG(3)^(3/2) / pi / sqrt((k(3)-1))/stabilise;

% eg basX*(VX2.*randn(size(VX2))) should have a mean RMS amplitude of 1;
% some type of Hanning like wrapper may be used in the future
VX2 = ones(k(1),1) * sqrt(VG(1)/k(1)) / stabilise;
VY2 = ones(k(2),1) * sqrt(VG(2)/k(2)) / stabilise;
VZ2 = ones(k(3),1) * sqrt(VG(3)/k(3)) / stabilise;

% multiply together
IC0 = [kron(VZ2,kron(VY2,VX1)); kron(VZ2,kron(VY1,VX2)); kron(VZ1,kron(VY2,VX2))].^(-2);
%IC0 = diag(IC0.^(-2));
%IC0 = [[IC0 zeros(size(IC0,1),size(VG,2))] ; zeros(size(VG,2),size(IC0,1)+size(VG,2))];
IC0 = sparse(1:length(IC0),1:length(IC0),IC0,length(IC0)+size(VG,2),length(IC0)+size(VG,2));

s1 = 3*prod(k);
s2 = s1 + size(VG,2);

T = zeros(s2,1);
T(s1+(1:size(VG,2))) = ones(size(VG,2),1);
for iter=1:param(4)
	fprintf('iteration # %d: ', iter);
	[Alpha,Beta,Var] = spm_brainwarp(VG,VF,Affine,basX,basY,basZ,T,fwhm);
	if (iter > 0)
		T = T + (Alpha/Var + IC0)\(Beta/Var - IC0*T);
	else
		T = T + Alpha\Beta;
	end
disp(Var);
end


remainder = T((1:size(VG,2)) + s1);
Transform = reshape(T(1:s1),prod(k),3)*stabilise.^3;
Dims = [VG(1:3,1)'; k];

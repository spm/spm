function [iy,y,z] = spm_pac_sample_deformation(z, fsubspace)
% FORMAT [iy,y] = sample_deformation([z0], [fname])
%   z0 - [M 1] Latent code. NaN values will be sampled from the
%        distribution {default: NaN}
%   fname - Path to the scaled subspace fail {'model/subspace_scaled.nii'}
%   iy - [Nx Ny Nz 3] Inverse deformation (used to forward deform meshes)
%   y  - [Nx Ny Nz 3] Forward deformation (used to forward deform volumes)
%   z  - [M 1] Sampled latent code

% The subspace that I gave you contains M = 100 principal components.
% This function assumes that the subspace contains eigenvectors scaled by
% the square root of their eigenvalues, so that sampling from the model can
% be done by sampling latent codes from the standard Normal distribution
% (mean:0, sd: 1)

%% Yael Balbastre 2024
% -------------------------------------------------------------------------
% Path to model file
if nargin < 2 || isempty(fsubspace)
    fsubspace = 'Templates/subspace_scaled.nii';
end

% -------------------------------------------------------------------------
% Load model
subspace = nifti(fsubspace);            % Memory map subspace in object U
subspace.dat.dim(4) = [];               % Remove (empty) 4th dimension
M   = size(subspace.dat,5);             % Number of principal components
dim = size(subspace.dat);               % Dimensions of the lattice
dim = dim(1:3);
vs = sqrt(sum(subspace.mat(1:3,1:3).^2));   % Voxel size
try
    % Parse regularisation parameters
    prm = strsplit(subspace.descrip, {' ' '(' ')'});
    prm = prm(2:6);
    prm = cellfun(@str2num, prm);
catch
    % Failed -> use shoot defaults
    dft = spm_shoot_defaults;
    prm = dft.rparam;
end

% -------------------------------------------------------------------------
% Sample latent code
if nargin < 1 || isempty(z)
    z = NaN(M,1);
end
z(~isfinite(z)) = randn(sum(~isfinite(z)),1);

% -------------------------------------------------------------------------
% Sample initial velocity
v = zeros([dim 3], 'single');
for m=1:M
    v = v + z(m) * single(subspace.dat(:,:,:,:,m));
end

% -------------------------------------------------------------------------
% Shoot deformations
if nargout == 1
    iy = spm_shoot3d(v, [vs prm], Inf);
else
    [iy,~,~,y] = spm_shoot3d(v, [vs prm], Inf);
end
function [iy,y,z] = spm_shp_sample_deformation(z,U,v0,z0)
% FORMAT [iy,y,z] = spm_shp_sample_deformation(z,U,v0,z0)
%
% z  - (M x 1) Latent code [NaN] 
%              NaN values will be sampled from the distribution 
% U  - Path to the scaled subspace [spm('Dir')/tpl/shp/subspace_scaled.nii]
% v0 - (Nx x Ny x Nz x 3) Original velocity field [zero]
% z0 - (Nx x Ny x Nz x 3) Original latent code    [zero]
% iy - (Nx x Ny x Nz x 3) Inverse deformation (used to deform meshes)
% y  - (Nx x Ny x Nz x 3) Forward deformation (used to deform volumes)
%
% This function:
% 1. Samples z values from the standard distribution (if nonfinite input)
% 2. Generates the velocity field v = v0 + U * (z - z0)
% 3. Exponentiates the forward and inverse deformation fields iy and y
%__________________________________________________________________________

% Yael Balbastre
% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging

% Note from email:
% The subspace that I gave you contains M = 100 principal components.
% This function assumes that the subspace contains eigenvectors scaled by
% the square root of their eigenvalues, so that sampling from the model can
% be done by sampling latent codes from the standard Normal distribution
% (mean:0, sd: 1)


% -------------------------------------------------------------------------
% Path to model file
if nargin < 1, z  = []; end
if nargin < 2, U  = []; end
if nargin < 3, v0 = []; end
if nargin < 4, z0 = []; end
if isempty(U) || (ischar(U) && ~exist(U, 'file'))
    U = spm_shp_get_model('subspace_scaled'); 
end


% -------------------------------------------------------------------------
% Load model
if ischar(U), U = nifti(U); end  % Memory map subspace in object U
U.dat.dim(4) = [];               % Remove (empty) 4th dimension
M   = size(U.dat,5);             % Number of principal components
dim = size(U.dat);               % Dimensions of the lattice
dim = dim(1:3);
vx  = sqrt(sum(U.mat(1:3,1:3).^2));   % Voxel size
try
    % Parse regularisation parameters
    prm = strsplit(U.descrip, {' ' '(' ')'});
    prm = prm(2:6);
    prm = cellfun(@str2num, prm);
catch
    % Failed -> use shoot defaults
    dft = spm_shoot_defaults;
    prm = dft.rparam;
end

% -------------------------------------------------------------------------
% Load velocity
if ischar(v0),       v0 = nifti(v0);                            end
if isa(v0, 'nifti'), v0 = v0.dat();                             end
if ~isempty(v0),     v0 = reshape(single(v0), [dim 3]);
else,                v0 = zeros([dim 3], 'single');             end
 
% -------------------------------------------------------------------------
% Sample latent code
if isempty(z0), z0 = zeros(M,1); end
if isempty(z),  z  = NaN(M,1);   end
z(~isfinite(z)) = randn(sum(~isfinite(z)),1);

% -------------------------------------------------------------------------
% Sample initial velocity
v = v0;
for m=1:M
    v = v + (z(m) - z0(m)) * single(U.dat(:,:,:,:,m));
end

% -------------------------------------------------------------------------
% Shoot deformations
if nargout == 1
    iy = spm_shoot3d(-v, [vx prm], Inf);
else
    [iy,~,~,y] = spm_shoot3d(-v, [vx prm], Inf);
end
function z = spm_shp_project_velocity(v, fmodel, fsubspace)
% Project a velocity (= a 4D volume) onto a subspace (= a 5D volume) to
% compute its latent code (a 1D vector).
%
% FORMAT z = spm_shp_project_velocity(v, [fmodel], [fsubspace])
%
% v         - (Nx x Ny x Nz x 3) Initial velocity
% fmodel    - Path to the model parameters
%             [spm('Dir')/tpl/shp/model_variables.mat]
% fsubspace - Path to the scaled subspace 
%             [spm('Dir')/tpl/shp/subspace_scaled.nii]
% z         - (M x 1) Latent code
%__________________________________________________________________________

% Yael Balbastre
% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Path to model file
if nargin < 2 || isempty(fmodel)
    fmodel = spm_shp_get_model('model_variables');
end
if nargin < 3 || isempty(fsubspace)
    fsubspace = spm_shp_get_model('subspace_scaled');
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
load(fmodel, 'A', 'Az', 'lam')
% A   = prior precision matrix over latent codes
% Az  = posterior precision matrix over latent codes
% lam = residual precision
sd = sqrt(1./diag(A));

% -------------------------------------------------------------------------
% Compute posterior latent code > E[z] = lam * inv(Az) * U' * L * v
    
% Compute L * v
Lv = spm_diffeo('vel2mom', single(v), double([vs prm]));

% Compute U' * (L*v)
z = zeros(M,1);
for m=1:M
    Um   = subspace.dat(:,:,:,:,m) / sd(m); % unscale PC
    z(m) = Um(:)' * Lv(:);
end
clear Lv Um

% Compute lam * inv(Az) * (U'*L*v)
z = lam*(Az\z);
z = z ./ sd; % scale latent
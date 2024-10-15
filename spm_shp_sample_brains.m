function [z,z0,outmesh] = spm_shp_sample_brains(varargin)
% FORMAT  [z,z0,mesh] = spm_shp_sample_brains(mesh, K, ...)
%
% Positional
% ----------
% mesh   - (N x 1) Mesh(es) to deform (gifti objects or paths).
% K      - Number of brains to make [default: 1]
%
% Keywords
% --------
% pc     - Vector of indices principal components to use [default: all]
% span   - Latent bound(s) [default: [0 ,3] ]
% fout   - Folder where to write generated gifti files [default: '.']
% fshp   - Folder where the shape model is stored.
% suffix - Output file suffix [default: 0]
% v0     - Subject's initial velocity [default: 0]
% y0     - Subject's transform [default: recompute]
% z0     - Subject's latent code [default: recompute]
% r2n    - Subject's import to native transform [default: identity]
% can    - If true:  center samples about canonical brain (z=0)
%          If false: center samples about subject's brain (z=z0)
%
% Returns
% -------
% z       - (M x K) Sampled latent codes
% z0      - (M x 1) Subject's latent code
% outmesh - (K x N) Deformed mesh(es) (gifti or paths)
%
% Some PC indices can be negative meaning their pca components (2) will be
% linearly shifted in opposite direction, i.e. [1 2] and [1 -2] are two
% different trajectories where 2nd component moves in opposite direction.
%______________________________  ____________________________________________

% Jose David Lopez, Yael Balbastre, Gareth Barnes
% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Input parser
p = inputParser;
p.addRequired('mesh');
p.addOptional('K',       1);
p.addParameter('pc',     [],   @isvector);
p.addParameter('span',   [0 3],    @isvector);
p.addParameter('fout',   [],   @ischar);
p.addParameter('fshp',   [],   @ischar);
p.addParameter('suffix', '0',  @ischar);
p.addParameter('can',    true, @(X) islogical(X) && isscalar(X));
p.addParameter('r2n',    []);
p.addParameter('v0',     []);
p.addParameter('y0',     []);
p.addParameter('z0',     []);
p.addParameter('zlimit',[]);
p.parse(varargin{:});

mesh    = p.Results.mesh;
K       = p.Results.K;
pc      = p.Results.pc;
span    = p.Results.span;
fout    = p.Results.fout;
fshp    = p.Results.fshp;
suffix  = p.Results.suffix;
v0      = p.Results.v0;
y0      = p.Results.y0;
z0      = p.Results.z0;
r2n     = p.Results.r2n;
can     = p.Results.can;
zlimit  = p.Results.zlimit;
% -------------------------------------------------------------------------
% Prepare inputs

% Load meshes as gifti
mesh0 = mesh;
if ~strcmpi(class(mesh), 'gifti')
    if ~iscell(mesh)
        mesh  = cellstr(mesh);
        mesh0 = mesh;
    end
    mesh = cellfun(@gifti, mesh, 'UniformOutput', false);
end

% Get model files
fmodel    = spm_shp_get_model('model_variables', fshp);
fsubspace = spm_shp_get_model('subspace_scaled', fshp);
subspace  = nifti(fsubspace);
mat0      = subspace.mat;
M         = subspace.dat.dim(6);

% Get initial velocity
if ~isempty(v0)
    % Memory map file
    if ischar(v0)
        v0 = nifti(v0);
    end
    % Load velocity
    if isa(v0, 'nifti')
        dim = v0.dat.dim(1:3);
        v0 = reshape(single(v0.dat()), [dim 3]);
    end
    % Project velocity
    if isempty(z0)
        z0 = spm_shp_project_velocity(v0, fmodel, fsubspace);
    end
    % Shoot velocity
    if isempty(y0)
        vx  = sqrt(sum(subspace.mat(1:3,1:3).^2));   % Voxel size
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
        [~,~,~,y0] = spm_shoot3d(-v0, [vx prm], Inf);
        y0 = spm_shp_warps('compose', mat0, y0);
    end
end

% Pre-transform model orientation matrix
if ~isempty(r2n)
    if ischar(r2n)
        load(r2n, 'r2n');
    end
else
    r2n = eye(4);
end
mat0 = r2n * mat0;

% Load subject-to-model transform + pre-transform mesh
if ~isempty(y0)
    if ischar(y0), y0 = nifti(y0); end
    if isa(y0, 'nifti')
        dim = y0.dat.dim(1:3);
        y0 = reshape(single(y0.dat()), [dim 3]);
    end
    y0 = spm_shp_warps('compose', r2n, y0);
    for n=1:numel(mesh)
        mesh{n} = spm_shp_transform_mesh(mesh{n}, y0, mat0);
    end
else
    for n=1:numel(mesh)
        mesh{n} = spm_shp_transform_mesh(mesh{n}, inv(r2n));
    end
end


% Load subject's code
if ischar(z0), load(z0, 'z'); z0 = z; end

% Create output directory
fsamp = fullfile(fout, 'Cerebros');
mkdir(fsamp);

% -------------------------------------------------------------------------
% Define sample space
% - Linearly spaced vector between [-sp sp] around z (subject) or zero (canonical)
% - Gaussian random distribution with mean over z (subject) or zero (canonical)

pc    = pc(:);      % ensure column vector
absPC = abs(pc);    % some components have negative values indicated we want to shift in opposite direction
sgnPC = sign(pc);   % + or -1 depending on how to shift
dz     = linspace(span(1),span(2),K);

dz = sgnPC .* dz;  % P x K




if ~can,
    z        = repmat(z0,1,K); %% make z the same as original brain by default
else
    z=zeros(M,K); %% poplation brain by default
end;

flip=0;



for j=1:length(absPC),
    for k=1:K,

        if (abs(z(absPC(j),k)+dz(j,k))>abs(zlimit)),
            dz(j,k)=-dz(j,k); 
            flip=flip+1;
        end; % if
        z(absPC(j),k)=z(absPC(j),k)+dz(j,k);
    end; % for k
end; % for j
fprintf('\n Flipped %3.2f percent of components to keep within abs(z)<%3.2f \n',100*flip/(K*length(absPC)),abs(zlimit));
% -------------------------------------------------------------------------
% Loop about samples
fprintf('Sample brain ');
outmesh = {};
for k=1:K
    if k > 1
        fprintf(repmat('\b', [1 4+numel(num2str(k-1))]));
    end
    fprintf('%d ...', k);

    % ---------------------------------------------------------------------
    % Sample transformation [template -> random]

    % Shoot along principal component
    z1 = z(:,k);


    % Build transformation field
    iy = spm_shp_sample_deformation(z1, subspace, v0, z0); % vox-to-vox
    iy = spm_shp_warps('compose', mat0, iy);               % vox-to-mm

    % ---------------------------------------------------------------------
    % Deform mesh
    for n=1:numel(mesh)

        % -----------------------------------------------------------------
        % Apply deformation - output is in template (MNI) space
        randomcortex = spm_shp_transform_mesh(mesh{n}, iy, mat0);

        % -----------------------------------------------------------------
        % Write to disk
        if ischar(mesh0{n})
            basename_mesh = spm_file(mesh0{n}, 'basename');
            name_mesh_nat = sprintf('smp_shp_%03d_%s%s.gii',k,basename_mesh,suffix);
            path_mesh_nat = fullfile(fsamp, name_mesh_nat);
            outmesh{k,n}  = path_mesh_nat;
            save(randomcortex, outmesh{k,n});
        else
            outmesh{k,n} = randomcortex;
        end
    end % for j
end % for n
fprintf('\n');



function [z0,zsubj,nativemeshnames] = spm_pca_sample_brains(Nb, PC, sp, output_folder, template_folder,template_meshfiles,randseed)
% FORMAT  spm_pca_sample_brains(Nb, PC, sp, output_folder, template_folder,template_meshfiles,randseed)
%   Nb - Number of brains to make
%   PC- vector of indices principal components to use eg [8:100]
% some PC indices can be negative meaning their pca components (2) will be linearly shifted in opposite direction
% i.e. [1 2] and [1 -2] are two different trajectories where 2nd component moves in opposite direction
%   sp- span in terms of z eg 3 gives range z=-3 -> +3
%   output_folder - Folder where to write generated gifti files {'.'}
%   template_folder- folder where the template brain principalc components
%   are stored 
% template_mesh_files- optional- defaults to cortex, but could also supply
% inner-skull mesh to be deformed etc
% rand_seed- optional unique file identifier (eg for different distortion
% trajectories)
%% Jose David Lopez, Yael Balbastre, Gareth Barnes 2024


%% Load data

% -------------------------------------------------------------------------
% Specify model files
if nargin<6,
    template_meshfiles=[];
end;
if nargin<7,
    randseed=[];
end;
if isempty(randseed),
    randseed=0;
end;

if isempty(template_meshfiles),
    template_meshfiles{1}=[output_folder 'mesh_cortex_template.gii'];
end;

if ~isvector(PC),
    error('PCin should be a vector (possible values 1:100)');
end;
fsubspace = [template_folder filesep 'subspace_scaled.nii'];
faffine   = [output_folder filesep 'affine.mat'];
mkdir([output_folder filesep 'Cerebros']);

% Subject components
fcode = [output_folder filesep 'latent_code.mat'];
load(fcode, 'z');						% JD: All PCs come from Subject's brain
zsubj=z; %% return original brain shape


% -------------------------------------------------------------------------
% Read model files
for f=1:numel(template_meshfiles), %% maybe more than 1 mesh eg cortex+ inner skull
    cortex{f}   = gifti([output_folder filesep template_meshfiles{f}]);
end;
subspace = nifti(fsubspace);
load(faffine, 'tpl2native');

%% Define way of randomisation
% - Creating a linearly spaced vector between [-sp sp] around z (subject) or zero (canonical)
% - Gaussian random distribution with mean over z (subject) or zero (canonical)

Npc = length(PC);
absPC=abs(PC); % some components have negative values indicated we want to shift in opposite direction
sgnPC=sign(PC); %% + or -1 depending on how to shift
if Npc == 1
    warning('moving one component systematically ')
    z0 = linspace(-sp,sp,Nb);				% Canonical
else
    fprintf('Linear spacing of lots of components..')
    lz = linspace(-sp,sp,Nb);				% Canonical
    z0=ones(Npc,Nb).*lz;
    z0=z0.*repmat(sgnPC',1,Nb);
end;

% -------------------------------------------------------------------------
% Loop about samples
fprintf('Sample brain ');
for n=1:Nb
    if n > 1
        fprintf(repmat('\b', [1 4+numel(num2str(n-1))]));
    end
    fprintf('%d ...', n);

    % ---------------------------------------------------------------------
    % Sample transformation [template -> random]

    % Shoot along principal component
    z = zeros(subspace.dat.dim(6),1);
    z(absPC) = z0(:,n);


    iy = spm_pca_sample_deformation(z, fsubspace);          % voxel-to-voxel
    iy = spm_pca_warps('compose', subspace.mat, iy);        % voxel-to-mm

    % ---------------------------------------------------------------------
    % Deform mesh
    for j=1:numel(template_meshfiles),
        [a1,b1,c1]=fileparts(template_meshfiles{j});
        namestr=[b1(1:findstr(b1,'template')-1) 'native']; %% output will be in native (MRI) space
        ncortex=cortex{j};
        randomcortex = spm_pca_transform_mesh(ncortex, iy, subspace.mat);

        % ---------------------------------------------------------------------
        % Send to Gareth space
        randomcortex = spm_pca_transform_mesh(randomcortex, tpl2native);

        % ---------------------------------------------------------------------
        % Write on disk
        nativemeshnames{n,j}=[output_folder filesep 'Cerebros\' sprintf('smp_%d_%s%d.gii',n,namestr,randseed)];
        save(randomcortex,nativemeshnames{n,j});
    end; % for j
end; % for n
fprintf('\n');



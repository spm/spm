function [G] = estimate_greens_mmclab(F, P)  
% Estimate Green's function of the photon fluence by simulating photon
% migration through the head and brain using the MMCLAB software
%
% This function provides input parameters to perform the Monte Carlo
% simulation. Output of MMCLAB, Green's function, is then saved as data
% format to be used in DCM-fNIRS analysis. 
% 
% Following software and Brain atlas are required: 
% MMCLAB: http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Doc/MMCLAB
% iso2mesh: http://iso2mesh.sourceforge.net/cgi-bin/index.cgi
% Brain atlas mesh: http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh
% Collin 27 average brain: http://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27
%
% FORMAT [G] = estimate_greens_mmclab(F, P) 
%
% G                 Green's function of the photon fluence 
%
% F                  names of files required to use MMClab software
% P                  optical parameters for Monte Carlo simulation
%
%--------------------------------------------------------------------------
% F.mmc          name of directory of MMCLAB (eg, /mmclab)
% F.isomesh     name of directory of iso2mesh (eg, /iso2mesh)
% F.mesh         file name of brain atlas mesh (eg, MMC_Collins_Atlas_Mesh_Version_2L.mat)
% F.atlas         file name of Collin 27 brain (eg, colin27_t1_tal_lin.nii) 
% F.sdpos        file name of MNI locations of optical source and detectors 
%                     If F is not specified, files are selected using GUI. 
%--------------------------------------------------------------------------
% P                  optical parameters for Monte Carlo simulation
%                     if this is not specified, optical parameters for 750
%                     nm and 850 nm are used. 
%--------------------------------------------------------------------------
% G.s - estimated Green's function from sensor (light emitter) positions
% into source positions [# sensor x # voxels x # wavelengths] 
% G.d - estimated Green's function from sensor (light detector) positions
% into source positions [# sensor x # voxels x # wavelengths] 
% G.xyz - MNI locations [3 x # voxels] 
% G.elem - tissue types of voxels [3 x # voxels] 
% 1-scalp, 2-CSF, 3-gray matter, 4-white matter 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: estimate_greens_mmclab.m 6754 2016-03-25 06:44:58Z will $

%-Load data and add folders of softwares to search path 
%--------------------------------------------------------------------------
if ~exist('F') 
    [F.mmc, sts(1)] = spm_select(1, 'dir', 'Select a directory of MMCLAB'); 
    [F.isomesh, sts(2)] = spm_select(1, 'dir', 'Select a directory of iso2mesh'); 
    [F.mesh, sts(3)] = spm_select(1, 'mat', 'Select a file of brain atlas mesh'); 
    [F.atlas, sts(4)] = spm_select(1, 'image', 'Select a file of Collin brain'); 
    [F.sdpos, sts(5)] = spm_select(1, 'mat', 'Select a file of MNI locations of optodes'); 
    
    if ~all(sts), 
        Gs = []; Gd = []; 
        return; 
    end
end

addpath(F.mmc); addpath(F.isomesh);
load(F.mesh); load(F.sdpos); 
vol = spm_vol(F.atlas);


%-Specify optical properties for Monte Carlo simulation 
%--------------------------------------------------------------------------
if ~exist('P') % initialize parameters 
    % Ref: AT Eggebrecht et al, NeuroImage 61 (2012): 1120-1128
    % absorption coefficient [750 nm 850nm] 
    mua(1,1) = 0.0143; mua(1,2) = 0.0164; % for scalp and skull
    mua(2,1) = 0.004;  mua(2,2) = 0.004; % CSF
    mua(3,1) = 0.018; mua(3,2) = 0.0192; % Gray matter
    mua(4,1) = 0.0167;  mua(4,2) = 0.0208; % White matter
    P.mua = mua; 
    
    % scattering coefficient [750 nm 850 nm] 
    mus(1,1) = 0.84; mus(1,2) = 0.74; 
    mus(2,1) = 0.3; mus(2,2) = 0.3;
    mus(3,1) = 0.8359; mus(3,2) = 0.6726;
    mus(4,1) = 1.1908; mus(4,2) = 1.0107;
    P.mus = mus; 
    
    % anisotropic factor 
    P.g = [0.89; 0.89; 0.89; 0.84]; 
    
    % refractory index
    P.n = ones(4,1) .* 1.37;
end
    
% reduced scattering coefficients 
rmus = P.mus ./ ((1-P.g)*ones(1,2));

prop(:,:,1) = [0 0 1 1; P.mua(:,1) rmus(:,1) P.g P.n]; 
prop(:,:,2) = [0 0 1 1; P.mua(:,2) rmus(:,2) P.g P.n]; 


%-Calculate initial positions and directions of photons to be launched 
%--------------------------------------------------------------------------
% transform the MNI coordinates of optical sources and detectors to the
% voxel coordinates of Colins27 Average Brain. 
vox_s = round(inv(vol.mat) * [pos.s.xyzH; ones(1, pos.s.ns)])'; 
vox_d = round(inv(vol.mat) * [pos.d.xyzH; ones(1, pos.d.nd)])'; 
vox_s(:,4) = []; vox_d(:,4) = []; 

[srcpos_s, srcdir_s] = find_start_photonpos(node, elem, face, vox_s); 
[srcpos_d, srcdir_d] = find_start_photonpos(node, elem, face, vox_d); 


%-Simulate photon migration through head and brain using MMCLab
%--------------------------------------------------------------------------
cfg.nphoton = 1e8; % the number of photons for MC simulation

cfg.tstart=0; % starting time of the simulation [seconds]
cfg.tend=3e-9; % ending time of the simulation [seconds]
cfg.tstep=2e-10; % time-gate width of the simulation [seconds]
cfg.debuglevel='TP';
cfg.node = node;
cfg.elem = elem;

% launch photons on source positions
for i = 1:2 % use two wavelengths of lights
    cfg.prop = prop(:,:,i);
    for j = 1:pos.s.ns
        cfg.srcpos = srcpos_s(j,:); 
        cfg.srcdir = srcdir_s(j,:); 
        out = mmclab(cfg); 
        
        flx = out.data; 
        G.s(j,:,i) = sum(flx, 2) .* cfg.tstep; 
    end
end

% launch photons on detector positions 
for i = 1:2 
    cfg.prop = prop(:,:,i); 
    for j = 1:pos.d.nd 
        cfg.srcpos = srcpos_d(j,:); 
        cfg.srcdir = srcdir_d(j,:); 
        out = mmclab(cfg); 
        
        flx = out.data .* 100; % unit: 1/(s*cm^2)
        G.d(j,:,i) = sum(flx, 2) .* cfg.tstep;
    end
end

G.xyz = vol.mat * [node'; ones(1, size(node, 1))]; G.xyz(4,:) = []; 
G.elem = elem; 


function [srcpos, srcdir] = find_start_photonpos(node, elem, face, vox_pos)
% Calculate initial positions and directions of photons to be launched 
% Note: dim(vox_pos) = [# voxels x 3] 

% calculate the centroid of the tetrahedral volume element
indx = elem(:, 1:4); 
xyz_e = node(indx(:),:); 
xyz_e = reshape(xyz_e, [size(elem, 1) 4 3]); 
cent_e = squeeze(mean(xyz_e, 2)); 

% calculate the centroid for the tetrahedral volume face
indx = face(:, 1:3); 
xyz_f = node(indx(:),:); 
xyz_f = reshape(xyz_f, [size(face, 1), 3 3]); 
cent = squeeze(mean(xyz_f, 2)); 

% find an initial position where photon is launched, 
% nearest centroid of tetrahedral volume element is selected. 
nvox = size(vox_pos, 1); 
for i = 1:nvox 
    dist = cent_e - ones(size(elem, 1), 1) * vox_pos(i,:);
    dist = sum(dist.^2, 2);
    indx = find(dist == min(dist));
    srcpos(i,:) = cent_e(indx,:);
end

% find an inital direction of photon movement 
% normal vector to tetrahedral volume face is calculated. 
AB=node(face(:,2),1:3)-node(face(:,1),1:3);
AC=node(face(:,3),1:3)-node(face(:,1),1:3);

N=cross(AB',AC')';
snorm = N./repmat(sqrt(sum(N.*N,2)),1,3);

for i = 1:nvox
    dist = cent - ones(size(face,1), 1) * srcpos(i,:);
    dist = sum(dist.^2, 2);
    indx = find(dist == min(dist));
    srcdir(i,:) = snorm(indx,:);
end

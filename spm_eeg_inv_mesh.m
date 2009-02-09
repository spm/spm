function mesh = spm_eeg_inv_mesh(sMRI, Msize)
% Apply the inverse spatial deformation to the template mesh
% to obtain the individual cortical mesh
% save the individual .mat tesselation of the chosen size
%
% FORMAT [fid, mesh] = spm_eeg_inv_meshing(filename, Msize)
% Input:
% sMRI - name of the sMRI file
% Msize - size of the mesh (1-3)
% Output:
% fid    - fiducials (head surface + points inverse normalized from the template)
% mesh   - inverse - normalized canonical mesh
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_mesh.m 2720 2009-02-09 19:50:46Z vladimir $

mesh = [];

% SPM directory of canonical anatomy
%--------------------------------------------------------------------------
Cdir     = fullfile(spm('dir'), 'canonical');

if nargin == 0 || isempty(sMRI)
    % Use the template
    %--------------------------------------------------------------------------
    mesh.template   = 1;
    mesh.Affine     = eye(4);
    mesh.sMRI       = fullfile(Cdir, 'single_subj_T1.nii');
else
    mesh.template = 0;
    mesh.sMRI = sMRI;

    mesh  = spm_eeg_inv_spatnorm(mesh);
    % ====================================
end

[pth, nam] = spm_fileparts(mesh.sMRI);

if nargin <2
    Msize = 2;
end

mesh.Msize = Msize;

spm('Pointer','Watch');

% Canonical cortical mesh
%------------------------------------------------------------------
switch mesh.Msize
    case 1
        filename = fullfile(Cdir, 'single_subj_cortex_5124.surf.gii');
    case 2
        filename = fullfile(Cdir, 'single_subj_cortex_8196.surf.gii');
    case 3
        filename = fullfile(Cdir, 'single_subj_cortex_20484.surf.gii');
end

mesh.tess_mni  = export(gifti(filename), 'spm');

% Compute the cortical mesh from the template
%------------------------------------------------------------------
if ~mesh.template
    Tmesh   = spm_swarp(filename, mesh.def);
    suffind = max(strfind(filename, '_cortex_'));
    filename = fullfile(pth, [nam filename(suffind:end)]);
    save(gifti(Tmesh), filename);
end

mesh.tess_ctx = filename;

% Compute the scalp mesh from the template
%------------------------------------------------------------------
filename = fullfile(Cdir, 'single_subj_scalp_2562.surf.gii');

if ~mesh.template
    Tmesh   = spm_swarp(filename, mesh.def);
    suffind = max(strfind(filename, '_scalp_'));
    filename = fullfile(pth, [nam filename(suffind:end)]);
    save(gifti(Tmesh), filename);
end

mesh.tess_scalp = filename;


% Compute the outer skull mesh from the template
%------------------------------------------------------------------
filename = fullfile(Cdir, 'single_subj_oskull_2562.surf.gii');

if ~mesh.template
    Tmesh   = spm_swarp(filename, mesh.def);
    suffind = max(strfind(filename, '_oskull_'));
    filename = fullfile(pth, [nam filename(suffind:end)]);
    save(gifti(Tmesh), filename);
end

mesh.tess_oskull = filename;

% Compute the inner skull mesh from the template
%------------------------------------------------------------------
filename = fullfile(Cdir, 'single_subj_iskull_2562.surf.gii');

if ~mesh.template
    Tmesh   = spm_swarp(filename, mesh.def);
    suffind = max(strfind(filename, '_iskull_'));
    filename = fullfile(pth, [nam filename(suffind:end)]);
    save(gifti(Tmesh), filename);
end

mesh.tess_iskull = filename;

% datareg
%--------------------------------------------------------------------------
mesh.fid = export(gifti(mesh.tess_scalp), 'ft');
mesh.fid.unit = 'mm';
mesh.fid.fid = struct('pnt', [1 85  -41;     -83 -20 -65;    83 -20 -65;    -87 -11 -62;     87 -11 -62], ...
    'label',{{'nas';          'lpa';         'rpa';       'FIL_CTF_L';    'FIL_CTF_R'}});

if ~mesh.template
    fidpnt = mesh.fid.fid;
    fidpnt = export(spm_swarp(gifti(fidpnt), mesh.def), 'ft');
    mesh.fid.fid.pnt = fidpnt.pnt;
end


spm('Pointer','Arrow');

function D = spm_eeg_inv_meshing(varargin)

%==========================================================================
% Apply the inverse spatial deformation to the template mesh
% to obtain the individual cortical mesh
% save the individual .mat tesselation of the chosen size
%
% FORMAT D = spm_eeg_inv_meshing(D,val)
% Input:
% D		    - input data struct (optional)
% Output:
% D			- same data struct including the new files and parameters
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_meshing.m 716 2007-01-16 21:13:50Z karl $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% ensure requisite field are specified
%--------------------------------------------------------------------------
if ~isfield(D.inv{val}.mesh,'sMRI')
    D     = spm_eeg_inv_spatnorm(D);
end

% get cortical mesh size
%--------------------------------------------------------------------------
if ~isfield(D.inv{val}.mesh,'Msize')
    Msize = spm_input('Mesh size (vertices)','+1','3000|4000|5000|7200',[1 2 3 4]);
    D.inv{val}.mesh.Msize = Msize;
else
    Msize = D.inv{val}.mesh.Msize;
end

% get or set canonical flag
%--------------------------------------------------------------------------
if ~isfield(D.inv{val}.mesh,'canonical')
    canonical = 1;
    D.inv{val}.mesh.canonical = canonical;
else
    canonical = D.inv{val}.mesh.canonical;
end

% set template flag
%--------------------------------------------------------------------------
D.inv{val}.mesh.template = 0;

% Compute the masks
%==========================================================================
D  = spm_eeg_inv_getmasks(D);

% Compute the skull, cortex and scalp meshes de nove, from masks
%==========================================================================
if ~canonical
    
    D  = spm_eeg_inv_getmeshes(D);

% or deform canonical meshes
%--------------------------------------------------------------------------
else

    % Compute the cortex mesh from the template
    %----------------------------------------------------------------------
    switch Msize
        case 1
            Tmesh = load('wmeshTemplate_3004d.mat');
        case 2
            Tmesh = load('wmeshTemplate_4004d.mat');
        case 3
            Tmesh = load('wmeshTemplate_5004d.mat');
        case 4
            Tmesh = load('wmeshTemplate_7204d.mat');
    end

    % Canonical cortical mesh
    %----------------------------------------------------------------------
    D.inv{val}.mesh.tess_mni.vert    = Tmesh.vert;
    D.inv{val}.mesh.tess_mni.face    = uint16(Tmesh.face);
    
    % Compute the cortical mesh from the template
    %----------------------------------------------------------------------
    Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{val}.mesh.def);
    D.inv{val}.mesh.tess_ctx.vert    = Tmesh.vert;
    D.inv{val}.mesh.tess_ctx.face    = uint16(Tmesh.face);

    % Compute the scalp mesh from the template
    %----------------------------------------------------------------------
    Tmesh      = load('wmeshTemplate_scalp.mat');
    Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{val}.mesh.def);
    D.inv{val}.mesh.tess_scalp.vert  = Tmesh.vert;
    D.inv{val}.mesh.tess_scalp.face  = uint16(Tmesh.face);

    % Compute the skull mesh from the template
    %----------------------------------------------------------------------
    Tmesh      = load('wmeshTemplate_skull.mat');
    Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{val}.mesh.def);
    D.inv{val}.mesh.tess_iskull.vert = Tmesh.vert;
    D.inv{val}.mesh.tess_iskull.face = uint16(Tmesh.face);

end
spm('Pointer','Arrow');

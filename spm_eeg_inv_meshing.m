function D = spm_eeg_inv_meshing(varargin)

%==========================================================================
% Apply the inverse spatial deformation to the template mesh
% to obtain the individual cortical mesh
% save the individual .mat tesselation of the chosen size
%
% FORMAT D = spm_eeg_inv_meshing(D,val)
% Input:
% D		   - input data struct (optional)
% Output:
% D		   - same data struct including the new files and parameters
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_meshing.m 934 2007-10-05 12:36:29Z karl $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% check template flag
%--------------------------------------------------------------------------
try
    template = D.inv{val}.mesh.template;
catch
    D.inv{val}.mesh.template = 0;
end

% get sMRI file name (select none if none present, to use template mesh)
%--------------------------------------------------------------------------
if ~D.inv{val}.mesh.template & ~isfield(D.inv{val}.mesh,'sMRI')
    D.inv{val}.mesh.sMRI = spm_select([0 1],'image','Select subject''s structural MRI (Press Done if none)');
end

if isempty(D.inv{val}.mesh.sMRI)
    disp('No structural MRI selected; will use template mesh');
    D.inv{val}.mesh.template = 1;
end

% set canonical flag (if not specified, determined by modality)
%--------------------------------------------------------------------------
if ~D.inv{val}.mesh.template
    try
        canonical = D.inv{val}.mesh.canonical;
    catch

        if strcmp(D.modality,'MEG')		    % No ECD yet for MEG!
            D.inv{val}.method = 'Imaging';
        end

        if ~isfield(D.inv{val},'method')
            D.inv{val}.method = questdlg('recontruction','Please select','Imaging','ECD','Imaging');
        end

        if strcmp(D.inv{val}.method,'ECD')	% Use subject-specific mesh for ECD
            D.inv{val}.mesh.canonical = 0;
        else
            D.inv{val}.mesh.canonical = 1;	% Use canonical mesh for Imaging
        end
    end
else
    D.inv{val}.mesh.canonical = 0;
end

% get cortical mesh size
%--------------------------------------------------------------------------
if ~isfield(D.inv{val}.mesh,'Msize')
    Msize = spm_input('Mesh size (vertices)','+1','3000|4000|5000|7200',[1 2 3 4]);
    D.inv{val}.mesh.Msize = Msize;
else
    Msize = D.inv{val}.mesh.Msize;
end


if D.inv{val}.mesh.template
    D = spm_eeg_inv_template(D);
else

    % Segment and normalise structural
    %======================================================================
    D  = spm_eeg_inv_spatnorm(D);

    % Compute the masks
    %======================================================================
    D  = spm_eeg_inv_getmasks(D);

    % Compute the skull, cortex and scalp meshes de novo, from masks
    %======================================================================
    if ~D.inv{val}.mesh.canonical

        D  = spm_eeg_inv_getmeshes(D);

        % or deform canonical meshes
        %------------------------------------------------------------------
    else

        % Compute the cortex mesh from the template
        %------------------------------------------------------------------
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
        %------------------------------------------------------------------
        D.inv{val}.mesh.tess_mni.vert    = Tmesh.vert;
        D.inv{val}.mesh.tess_mni.face    = uint16(Tmesh.face);

        % Compute the cortical mesh from the template
        %------------------------------------------------------------------
        Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{val}.mesh.def);
        D.inv{val}.mesh.tess_ctx.vert    = Tmesh.vert;
        D.inv{val}.mesh.tess_ctx.face    = uint16(Tmesh.face);

        % Compute the scalp mesh from the template
        %------------------------------------------------------------------
        Tmesh      = load('wmeshTemplate_scalp.mat');
        Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{val}.mesh.def);
        D.inv{val}.mesh.tess_scalp.vert  = Tmesh.vert;
        D.inv{val}.mesh.tess_scalp.face  = uint16(Tmesh.face);

        % Compute the skull mesh from the template
        %------------------------------------------------------------------
        Tmesh      = load('wmeshTemplate_skull.mat');
        Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{val}.mesh.def);
        D.inv{val}.mesh.tess_iskull.vert = Tmesh.vert;
        D.inv{val}.mesh.tess_iskull.face = uint16(Tmesh.face);
    end
end
spm('Pointer','Arrow');

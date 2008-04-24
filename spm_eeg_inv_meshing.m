function D = spm_eeg_inv_meshing(varargin)
% Apply the inverse spatial deformation to the template mesh
% to obtain the individual cortical mesh
% save the individual .mat tesselation of the chosen size
%
% FORMAT D = spm_eeg_inv_meshing(D,ival)
% Input:
% D        - input data struct (optional)
% Output:
% D        - same data struct including the new files and parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_meshing.m 1477 2008-04-24 14:33:47Z christophe $

% Flags explanation:
% * template:
%   [1] - use the template mesh, useful when the subject's sMRI is not
%         available
%   [0] - use the subject's sMRI to estimate at least the main "brain
%         volumes", and maybe other meshes.
% * canonical:
%   [2] - use the canonical meshes for the cortex, i/o skull and scalp and 
%         bring it into the subject's own space, thanks to the 
%         normalisation parameters of his sMRI (if available).
%   [1] - use the canonical cortiex mesh and subject's specifi meshes for
%         the inner skull and scalp.
%   [0] - extract the subject's meshes from the binary volumes estimated
%         from his sMRI.
%         Beware that the cortex mesh will be a very rough approximation !
%         A mesh estimated by a more sophisticated toolbox (eg. BrainVisa)
%         would provide a more accurate mesh.


% initialise
%--------------------------------------------------------------------------
[D,ival] = spm_eeg_inv_check(varargin{:});

% check template flag
%--------------------------------------------------------------------------
try
    template = D.inv{ival}.mesh.template;
catch
    D.inv{ival}.mesh.template = 0;
end

% get sMRI file name (select none if none present, to use template mesh)
%--------------------------------------------------------------------------
if ~D.inv{ival}.mesh.template & ~isfield(D.inv{ival}.mesh,'sMRI')
    D.inv{ival}.mesh.sMRI = spm_select([0 1],'image', ...
                'Select subject''s structural MRI (Press Done if none)');
end

if isempty(D.inv{ival}.mesh.sMRI)
    disp('No structural MRI selected; will use template mesh');
    D.inv{ival}.mesh.template = 1;
end

% set canonical flag (if not specified, determined by modality)
%--------------------------------------------------------------------------
if ~D.inv{ival}.mesh.template
    try
        canonical = D.inv{ival}.mesh.canonical;
    catch

        if strcmp(D.modality,'MEG')         % No ECD yet for MEG!
            D.inv{ival}.method = 'Imaging';
        end

        if ~isfield(D.inv{ival},'method')
            D.inv{ival}.method = questdlg('recontruction','Please select', ...
                                                'Imaging','ECD','Imaging');
        end

        if strcmp(D.inv{ival}.method,'ECD')  % Use subject-specific mesh for ECD
            D.inv{ival}.mesh.canonical = 1;
            % We could allow to use template specific model for ECD by
            % allowing canonical = 2;
        else % Imaging
            if ~isfield(D.inv{ival}.mesh,'canonical')
                tmp = questdlg('scalp/skull mesh','Please select', ...
                                'canonical','subject','canonical');
            end
            if strcmp(tmp,'canonical')
                D.inv{ival}.mesh.canonical = 2;  % Use canonical mesh
            else
                D.inv{ival}.mesh.canonical = 1;
            end
        end
    end
else
    D.inv{ival}.mesh.canonical = 2;
end
% The case with canonical=0, i.e. use all meshes (scalp/iskull/cortex) 
% defined from the subject's sMRI, only occurs if canonical is defined by
% hand somewhere else.

% get cortical mesh size
%--------------------------------------------------------------------------
if ~isfield(D.inv{ival}.mesh,'Msize') ...
                && strcmp(D.inv{ival}.method,'Imaging')
    Msize = spm_input('Cortical mesh size (vert.)','+1', ...
                                          '3000|4000|5000|7200',[1 2 3 4]);
    D.inv{ival}.mesh.Msize = Msize;
else strcmp(D.inv{ival}.method,'ECD')
    D.inv{ival}.mesh.Msize = 1;
end
Msize = D.inv{ival}.mesh.Msize;

if D.inv{ival}.mesh.template
    [vol,fid,mesh] = spm_eeg_inv_template(Msize);
    D = fiducials(D,fid);
    D.inv{ival}.mesh = mesh;
    D.inv{ival}.vol = vol;
else
    mesh = D.inv{ival}.mesh;
    % Segment and normalise structural
    %======================================================================
    mesh  = spm_eeg_inv_spatnorm(mesh,ival);
    
    % Compute the masks
    %======================================================================
    mesh  = spm_eeg_inv_getmasks(mesh);
    D.inv{ival}.mesh = mesh;

    % Compute the skull, cortex and scalp meshes de novo, from masks
    %======================================================================
    if mesh.canonical<2
        
        switch mesh.canonical
            case 1
                [mesh,vol] = spm_eeg_inv_getmeshes(mesh,1:3);
            case 2
                [mesh,vol] = spm_eeg_inv_getmeshes(mesh,1:4);
        end

        % or deform canonical meshes
        %------------------------------------------------------------------
    else
        % Get mesh from template and warp it into subject's space
        %------------------------------------------------------------------
        tmp = load('standard_SPMvol.mat'); % rebuilt with SPM
%         tmp = load('standard_vol.mat');  % original from FT

        vol_mni = tmp.vol;
        vol = vol_mni;
        % deal with the meshes
        %------------------------------------------------------------------
        for ii=1:length(vol.bnd)
            vol.bnd(ii).pnt = spm_get_orig_coord(vol.bnd(ii).pnt,mesh.def);
        end
% 
%         % Compute the scalp mesh from the template
%         %------------------------------------------------------------------
%         Tmesh      = load('wmeshTemplate_scalp.mat');
%         Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{ival}.mesh.def);
%         D.inv{ival}.mesh.tess_scalp.vert  = Tmesh.vert;
%         D.inv{ival}.mesh.tess_scalp.face  = uint16(Tmesh.face);
%         % Compute the skull mesh from the template
%         %------------------------------------------------------------------
%         Tmesh      = load('wmeshTemplate_skull.mat');
%         Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{ival}.mesh.def);
%         D.inv{ival}.mesh.tess_iskull.vert = Tmesh.vert;
%         D.inv{ival}.mesh.tess_iskull.face = uint16(Tmesh.face);
    end
    D.inv{ival}.vol = vol;

    % Deal with the cortical mesh, canonical cortex mesh if 1 or 2
    if mesh.canonical
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
        mesh.tess_mni.vert    = Tmesh.vert;
        mesh.tess_mni.face    = uint16(Tmesh.face);

        % Compute the cortical mesh from the template
        %------------------------------------------------------------------
        Tmesh.vert = spm_get_orig_coord(Tmesh.vert,D.inv{ival}.mesh.def);
        mesh.tess_ctx.vert    = Tmesh.vert;
        mesh.tess_ctx.face    = uint16(Tmesh.face);
    end
    D.inv{ival}.mesh = mesh;
end
spm('Pointer','Arrow');

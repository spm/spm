function [D, L] = spm_opm_headmodel(S)
% Coregister FIL OPM data and option to set up forward model
% FORMAT D = spm_opm_create(S)
%   S               - input structure
% Fields of S:
%   S.D             - SPM MEEG Object                          - Default: REQUIRED 
%   S.coordsystem   - coordsystem.json file                    - Default: transform between sensor space and anatomy is identity
%   S.sMRI          - Filepath to  MRI file                    - Default: no Default
%   S.template      - Use SPM canonical template               - Default: 0
%   S.headhape      - .pos file for better template fit        - Default: no Default
%   S.cortex        - Custom cortical mesh                     - Default: Use inverse normalised cortical mesh
%   S.scalp         - Custom scalp mesh                        - Default: Use inverse normalised scalp mesh
%   S.oskull        - Custom outer skull mesh                  - Default: Use inverse normalised outer skull mesh
%   S.iskull        - Custom inner skull mesh                  - Default: Use inverse normalised inner skull mesh
%   S.voltype       - Volume conducter Model type              - Default: 'Single Shell'
%   S.meshres       - mesh resolution(1,2,3)                   - Default: 1
%   S.lead          - flag to compute lead field               - Default: 0
% Output:
%  D           - MEEG object (also written to disk)
%  L           - Lead field (also written on disk)
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


spm('FnBanner', mfilename);

if ~isfield(S,'D'); error('please specify a MEEG object!'); end
if ~isfield(S,'template');      S.template = 0;             end
if ~isfield(S,'sMRI');          S.sMRI = [];                end
if ~isfield(S,'coordsystem');   S.coordsystem = [];         end
if ~isfield(S,'headshape');     S.headshape = [];           end
if ~isfield(S,'meshres');       S.meshres = 1;              end
if ~isfield(S,'iskull');        S.iskull = [];              end
if ~isfield(S,'oskull');        S.oskull = [];              end
if ~isfield(S,'scalp');         S.scalp = [];               end
if ~isfield(S,'cortex');        S.cortex = [];              end
if ~isfield(S,'voltype');       S.voltype = 'Single Shell'; end
if ~isfield(S,'lead');          S.lead = 0;                 end

D = spm_eeg_load(S.D);

% check if template brain to be used
if S.template && isempty(S.sMRI)
    S.sMRI = 1;
elseif isnumeric(S.sMRI)
    if S.sMRI
        S.template = 1;
    end
end

% check whether coordsystem actually exists, as it might have been
% auto-imported from spm_opm_create
if ~isempty(S.coordsystem)
    if ~exist(S.coordsystem,'file')
        warning('coordsystem not found');
        S.coordsystem = [];
    end
end

% identify what kind of headmodel specification to do, this should result
% in a mutually exclusive job dictated
do.nothing              = 0;
do.fiducials            = 0;    % Just add fiducial field
do.identity             = 0;    % sensors and MRI are already in same space
do.coreg                = 0;

if isempty(S.coordsystem) && isempty(S.sMRI)
    do.nothing = 1;
elseif ~isempty(S.coordsystem) && isempty(S.sMRI)
    do.fiducials = 1;
elseif isempty(S.coordsystem) && ~isempty(S.sMRI)
    do.identity = 1;
elseif ~isempty(S.coordsystem) && ~isempty(S.sMRI)
    do.coreg = 1;
end

% deal with the nothing job
if do.nothing
    warning('No files for coregistration supplied, skipping\n')
end

%- MRI and meshes.
%--------------------------------------------------------
if do.identity || do.coreg
    %initially used inverse normalised meshes
    D = spm_eeg_inv_mesh_ui(D,1,S.sMRI,S.meshres);
    save(D);
    % then fill in custom meshes(if they exist)
    args = [];
    args.D = D;
    args.scalp = S.scalp;
    args.cortex = S.cortex;
    args.iskull = S.iskull;
    args.oskull = S.oskull;
    args.template = S.template;
    D = opm_customMeshes(args);
    save(D);
end

%- MEG Fiducials
%--------------------------------------------------------
fid = [];
if do.fiducials || do.coreg
    
    coord = spm_load(S.coordsystem);
    hnames = fieldnames(coord.HeadCoilCoordinates);
    fiMat = zeros(numel(hnames),3);
    for ii = 1:numel(hnames)
        fiMat(ii,:) = coord.HeadCoilCoordinates.(hnames{ii});
    end
    fid.fid.label = hnames;
    fid.fid.pnt = fiMat;
    % load in headshape if supplied
    if isempty(S.headshape)
        fid.pnt = [];
    else
        % SPMs native support for polhemus files was ditched a long
        % time ago, so using fieldTrip as a backend.
        shape = ft_read_headshape(S.headshape);
        fid.pnt = shape.pos;
    end
    
elseif do.identity
    warning('Assuming sensors and MRI are in same space\n')
    fid.fid.label = {'nas', 'lpa', 'rpa'}';
    fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:); % George O'Neill
    fid.pnt = [];
end

if ~isempty(fid)
    D = fiducials(D, fid);
    save(D);
end

%- MRI Landmark specification
%--------------------------------------------------------
if do.coreg
    M = [];
    % check if MRI landmarks are supplied in coordsystem file
    if isfield(coord,'AnatomicalLandmarkCoordinates')
        anames = fieldnames(coord.AnatomicalLandmarkCoordinates);
        % check all fields from hnames are in anames
        if ~isempty(setdiff(hnames,anames))
            error(['Some fidcuial coils are not mentioned in '...
                'the anatomical landmark description. Please '...
                'check your coordsystem.json file.']);
        end
        M.fid.label = hnames;
        M.fid.pnt = zeros(numel(hnames),3);
        for ii = 1:numel(hnames)
            M.fid.pnt(ii,:) = coord.AnatomicalLandmarkCoordinates.(hnames{ii});
        end
        M.pnt = [];
    else
        warning(['Anatomical landmarks not found in coordsystem.json, '...
            'will use derived fiducials from MRI itself']);
        
        % try and match up with some common names
        targets = {'nz','nas','nasion','lpa','rpa'};
        missing = setdiff(lower(hnames),targets);
        if ~isempty(missing)
            warning(['fiducial naming unclear, '...
                'assuming nas/lpa/rpa ordering']);
            M.fid.label     = {'nas', 'lpa', 'rpa'}';
            M.fid.pnt       = D.inv{1}.mesh.fid.fid.pnt(1:3,:);
            
        else
            for ii = 1:numel(hnames)
                switch find(strcmpi(targets,hnames{ii}))
                    case {1,2,3}
                        % nasion
                        M.fid.label{ii} = 'nas';
                        M.fid.pnt(ii,:) = D.inv{1}.mesh.fid.fid.pnt(1,:);
                    case 4
                        % lpa
                        M.fid.label{ii} = 'lpa';
                        M.fid.pnt(ii,:) = D.inv{1}.mesh.fid.fid.pnt(2,:);
                    case 5
                        % rpa
                        M.fid.label{ii} = 'rpa';
                        M.fid.pnt(ii,:) = D.inv{1}.mesh.fid.fid.pnt(3,:);
                end
            end
        end
        if isempty(S.headshape)
            M.pnt           = [];
        else
            M.pnt = D.inv{1}.mesh.fid.pnt;
        end
        
    end
    
elseif do.identity
    
    M.fid.label     = {'nas', 'lpa', 'rpa'}';
    M.fid.pnt       = D.inv{1}.mesh.fid.fid.pnt(1:3,:);
    if isempty(S.headshape)
        M.pnt           = [];
    else
        M.pnt = D.inv{1}.mesh.fid.pnt;
    end
    
end

%- Registration
%--------------------------------------------------------
if do.coreg || do.identity
    f = fiducials(D);
    if ~isempty(f.pnt)
        D = spm_eeg_inv_datareg_ui(D,1,f,M,1);
    else
        f.pnt = zeros(0,3);
        D = spm_eeg_inv_datareg_ui(D,1,f,M,0);
    end
end

save(D);

%- Foward model specification / calculation
%--------------------------------------------------------------------------
if do.coreg || do.identity
    
    D.inv{1}.forward.voltype = S.voltype;
    D = spm_eeg_inv_forward(D);
    if(S.lead)
        nverts = length(D.inv{1}.forward.mesh.vert);
        [L,D] = spm_eeg_lgainmat(D,1:nverts);
    end
    spm_eeg_inv_checkforward(D,1,1);
end
save(D);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Custom Meshes                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = opm_customMeshes(S)
% wrapper for adding custom meshes to MEG object
% FORMAT D = spm_opm_customMeshes(S)
%   S               - input structure
% Fields of S:
%   S.D             - Valid MEG object         - Default:
%   S.cortex        - Cortical mesh file       - Default: Use inverse normalised cortical mesh
%   S.scalp         - Scalp mesh file          - Default: Use inverse normalised scalp mesh
%   S.oskull        - Outer skull mesh file    - Default: Use inverse normalised outer skull mesh
%   S.iskull        - Inner skull mesh file    - Default: Use inverse normalised inner skull mesh
%   S.template      - is mesh in MNI space?    - Default: 0
% Output:
%  D           - MEEG object
%--------------------------------------------------------------------------

%- Default values & argument check
%--------------------------------------------------------------------------
if ~isfield(S, 'D'),           error('MEG object needs to be supplied'); end
if ~isfield(S, 'cortex'),      S.cortex = []; end
if ~isfield(S, 'scalp'),       S.scalp = []; end
if ~isfield(S, 'oskull'),      S.oskull = []; end
if ~isfield(S, 'iskull'),      S.iskull = []; end
if ~isfield(S, 'template'),    S.template = 0; end

D = S.D;
if ~isfield(D.inv{1}.mesh,'sMRI')
    error('MEG object needs to be contain inverse normalised meshes already')
end

%- add custom scalp and skull meshes if supplied
%--------------------------------------------------------------------------
if ~isempty(S.scalp)
    D.inv{1}.mesh.tess_scalp = S.scalp;
end

if ~isempty(S.oskull)
    D.inv{1}.mesh.tess_oskull = S.oskull;
end

if ~isempty(S.iskull)
    D.inv{1}.mesh.tess_iskull = S.iskull;
end

%- add custom cortex and replace MNI cortex with warped cortex
%--------------------------------------------------------------------------
if ~isempty(S.cortex)
    D.inv{1}.mesh.tess_ctx = S.cortex;
    if(S.template)
        D.inv{1}.mesh.tess_mni = S.cortex;
    else
        defs.comp{1}.inv.comp{1}.def = {D.inv{1}.mesh.def};
        defs.comp{1}.inv.space = {D.inv{1}.mesh.sMRI};
        defs.out{1}.surf.surface = {D.inv{1}.mesh.tess_ctx};
        defs.out{1}.surf.savedir.savesrc = 1;
        out = spm_deformations(defs);
        D.inv{1}.mesh.tess_mni  = export(gifti(out.surf{1}), 'spm');
    end
end
save(D);

end
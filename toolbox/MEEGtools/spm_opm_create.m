function [D,L] = spm_opm_create(S)
% Read magnetometer data and optionally set up forward model
% FORMAT D = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
% SENSOR LEVEL INFO
%   S.data          - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.channels      - channels.tsv file                        - Default: REQUIRED
%   S.fs            - Sampling frequency (Hz)                  - Default: REQUIRED if S.meg is empty
%   S.meg           - meg.json file                            - Default: REQUIRED if S.fs is empty
%   S.precision     - 'single' or 'double'                     - Default: 'single'
% SOURCE LEVEL INFO
%   S.coordsystem   - coordsystem.json file                    - Default: transform between sensor space and anatomy is identity
%   S.positions     - positions.tsv file                       - Default: no Default
%   S.sMRI          - Filepath to  MRI file                    - Default: no Default
%   S.template      - Use SPM canonical template               - Default: 0
%   S.headhape      - .pos file for better template fit        - Default:
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
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney
% $Id: spm_opm_create.m 8201 2021-12-22 14:11:18Z george $
spm('FnBanner', mfilename);

%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'voltype'),     S.voltype = 'Single Shell';  end
if ~isfield(S, 'meshres'),     S.meshres = 1;               end
if ~isfield(S, 'scalp'),       S.scalp = [];                end
if ~isfield(S, 'template'),    S.template = 0;              end
if ~isfield(S, 'cortex'),      S.cortex = [];               end
if ~isfield(S, 'iskull'),      S.iskull = [];               end
if ~isfield(S, 'oskull'),      S.oskull = [];               end
if ~isfield(S, 'fname'),       S.fname = 'sim_opm';         end
if ~isfield(S, 'path'),        S.path = [];                 end
if ~isfield(S, 'precision'),   S.precision = 'single';      end
if ~isfield(S, 'lead'),        S.lead = 0;                  end
if ~isfield(S, 'headshape');   S.headshape = [];            end


%- identify Binary File
%----------------------------------------------------------------------
try % work out if data is a matrix or a file
    [direc, dataFile] = fileparts(S.data);
    binData=1;
catch % if not readable check if it is numeric
    if ~isa(S.data,'numeric') % if not numeric throw error
        error('A valid dataest or file was not supplied')
    end
    binData=0;
    direc = pwd();
    dataFile=S.fname;
    if ~isfield(S, 'channels')
        error('A channels.tsv file must be supplied');
    end
end
%- identify potential BIDS Files
%----------------------------------------------------------------------
base = strsplit(dataFile,'meg');
chanFile= fullfile(direc,[base{1},'channels.tsv']);
megFile= fullfile(direc,[base{1},'meg.json']);
posFile= spm_select('FPList',direc,[base{1},'positions.tsv']);
coordFile= fullfile(direc,[base{1},'coordsystem.json']);

%- Check for channel Info
%--------------------------------------------------------------------------
try % to load a channels file
    channels = spm_load(S.channels);
catch
    try % to load a BIDS channel file
        channels = spm_load(chanFile);
    catch
        try  % use channel struct if supplied
            channels = S.channels;
        catch % create channel struct
            error('A valid channels.tsv file or struct was not found');
        end
    end
end

%- Check for MEG Info
%--------------------------------------------------------------------------
try % to load a meg file
    meg = spm_load(S.meg);
catch
    try % to load a BIDS meg file
        meg = spm_load(megFile);
    catch
        try % to use meg struct
            meg = S.meg;
        catch
            try % to use S.fs argument to get sampling frequency
                meg =[];
                meg.SamplingFrequency=S.fs;
            catch
                error('A meg.json file is required if S.fs is empty');
            end
        end
    end
end

%- Position File check
%----------------------------------------------------------------------
try % to load a channels file
    posOri = spm_load(S.positions);
    positions =1;
catch
    try % to load a BIDS channel file
        posOri = spm_load(posFile);
        positions =1;
    catch
        try % to assign a BIDS struct of positions
            if (isnumeric(S.positions))
                positions=1;
            else
                positions =0;
            end
        catch
            warning('No position information found')
            positions=0;
        end
    end
end

%- Forward model Check
%----------------------------------------------------------------------
subjectSource   = positions & isfield(S,'sMRI');
subjectNoStruct = positions & S.template == 1;

if subjectSource
    switch class(S.sMRI)
        case 'char'
            forward = 1;
            template = 0;
        case 'double'
            forward         = 1;
            template        = S.sMRI==1;
    end
elseif subjectNoStruct
    forward         = 1;
    template        = 1;
    S.sMRI          = 1;
else
    forward             = 0;
    template            = 0;
end

%- work out data size
%--------------------------------------------------------------------------
nChans = size(channels.name,1);

if(binData)
    fprops= dir(S.data);
    if(S.precision == 'single')
        bytesPerSample=4;
    else
        bytesPerSample=8;
    end
    nSamples = fprops.bytes/(nChans*bytesPerSample);
    nTrials = 1;
else
    nSamples=size(S.data,2);
    nTrials=size(S.data,3);
end

%- Check for a custom save path
%--------------------------------------------------------------------------
if ~isempty(S.path)
    if exist(S.path,'dir')
        direc = S.path;
    else
        error('specified output directory does not exist!')
    end
end % will use original direc variable otherwise

%- Create SPM object
%--------------------------------------------------------------------------

D = meeg(nChans,nSamples,nTrials);
D = fsample(D,meg.SamplingFrequency);
D = fname(D,[dataFile,'.mat']);
D = path(D,direc);
D = chanlabels(D,1:size(D,1),channels.name);
D = units(D,1:size(D,1),channels.units);
D = chantype(D,1:size(D,1),channels.type);

%- Overwrite and Save
%--------------------------------------------------------------------------
ma = fullfile(direc,[dataFile,'.mat']);
da = fullfile(direc,[dataFile,'.dat']);

ae = exist(fullfile(path(D),fname(D)),'file')==2;
if(ae)
    delete(ma);
    delete(da);
end
D.save();

%- reformat data according to channel info
%--------------------------------------------------------------------------
D= blank(D,[dataFile,'.dat']);
if binData
    
    maxMem= 100e6;
    samplesPerChunk= round(maxMem/(8*nChans));
    begs= 1:samplesPerChunk:nSamples;
    ends= begs+samplesPerChunk-1;
    
    if(ends(end)>nSamples)
        ends(end)=nSamples;
    end
    
    dat = fopen(S.data);
    for j=1:length(begs)
        asamplesPerChunk=length(begs(j):ends(j));
        inds= begs(j):ends(j);
        D(:,inds,1) = fread(dat,[nChans,asamplesPerChunk],S.precision,0,'b');
    end
    fclose(dat);
else
    % insert provided data
    D(1:nChans,1:nSamples,1:nTrials) = S.data;
    D.save();
end

%- Create Meshes
%--------------------------------------------------------------------------
if forward
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
    args.template = template;
    D = opm_customMeshes(args);
    save(D);
end

%-Place Sensors  in object
%--------------------------------------------------------------------------
if positions
    pos = [posOri.Px,posOri.Py,posOri.Pz];
    ori = [posOri.Ox,posOri.Oy,posOri.Oz];
    cl = posOri.name;
    
    [sel1 sel2] = match_str(cl,channels.name);
    
    grad= [];
    grad.label = cl;
    grad.coilpos = pos;
    grad.coilori = ori;
    grad.tra = eye(numel(grad.label));
    grad.chanunit = repmat({'T'}, numel(grad.label), 1);
    grad.chantype = lower({channels.type{sel2}}');
    grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'mm');
    D = sensors(D, 'MEG', grad);
    save(D);
    
    %- 2D view based on mean orientation of sensors
    n1=mean(grad.coilori); n1= n1./sqrt(dot(n1,n1));
    t1=cross(n1,[0 0 1]);
    t2=cross(t1,n1);
    pos2d =zeros(size(grad.coilpos,1),2);
    for i=1:size(grad.coilpos,1)
        pos2d(i,1)=dot(grad.coilpos(i,:),t1);
        pos2d(i,2)=dot(grad.coilpos(i,:),t2);
    end
    
    nMEG = length(indchantype(D,'MEG'));
    if nMEG~=size(pos2d,1)
        m1 = '2D positions could not be set as there are ';
        m2 =num2str(nMEG);
        m3 = ' channels but only ';
        m4 = num2str(size(pos2d,1));
        m5 =  ' channels with position information.';
        message = [m1,m2,m3,m4,m5];
        warning(message);
    else
        args=[];
        args.D=D;
        args.xy= pos2d';
        args.label=grad.label;
        args.task='setcoor2d';
        D=spm_eeg_prep(args);
        D.save;
    end
end


%- fiducial settings
%--------------------------------------------------------------------------
if subjectSource
    miMat = zeros(3,3);
    fiMat = zeros(3,3);
    fid=[];
    try % to read the coordsystem.json
        coord = spm_load(S.coordsystem);
        fiMat(1,:) = coord.HeadCoilCoordinates.coil1;
        fiMat(2,:) = coord.HeadCoilCoordinates.coil2;
        fiMat(3,:) = coord.HeadCoilCoordinates.coil3;
        miMat(1,:) = coord.AnatomicalLandmarkCoordinates.coil1;
        miMat(2,:) = coord.AnatomicalLandmarkCoordinates.coil2;
        miMat(3,:) = coord.AnatomicalLandmarkCoordinates.coil3;
        fid.fid.label = fieldnames(coord.HeadCoilCoordinates);
        fid.fid.pnt =fiMat;
        fid.pos= []; % headshape field that is left blank (GRB)
        M = fid;
        M.fid.pnt=miMat;
        M.pnt = D.inv{1}.mesh.fid.pnt;
    catch
        try % to find the BIDS coordsystem.json
            coord = spm_load(coordFile);
            fiMat(1,:) = coord.HeadCoilCoordinates.coil1;
            fiMat(2,:) = coord.HeadCoilCoordinates.coil2;
            fiMat(3,:) = coord.HeadCoilCoordinates.coil3;
            miMat(1,:) = coord.AnatomicalLandmarkCoordinates.coil1;
            miMat(2,:) = coord.AnatomicalLandmarkCoordinates.coil2;
            miMat(3,:) = coord.AnatomicalLandmarkCoordinates.coil3;
            fid.fid.label = fieldnames(coord.HeadCoilCoordinates);
            fid.fid.pnt =fiMat;
            fid.pos= []; % headshape field that is left blank (GRB)
            M = fid;
            M.fid.pnt=miMat;
            M.pnt = D.inv{1}.mesh.fid.pnt;
        catch % DEFAULT: transform between fiducials and anatomy is identity
            fid.fid.label = {'nas', 'lpa', 'rpa'}';
            fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:); % George O'Neill
            fid.pos= [];
            M = fid;
            M.pnt = D.inv{1}.mesh.fid.pnt; % George O'Neill
        end
    end
end

% If user wants to use the template brain, try to load fiducials from
% coordsystem.json
if subjectNoStruct
    try
        
        % check if there is a headshape file and load data/fids from there
        % otherwise fallback on coordsystem.json
        if ~isempty(S.headshape)
            % SPMs native support for polhemus files was ditched a long
            % time ago, so using fieldTrip as a backend.
            shape = ft_read_headshape(S.headshape);
            fid.pnt = shape.pos;
            
            % headshape mush be in the same space as sensors!
            % get the order into nas, lpa, rpa
            targetOrder = {'nas','lpa','rpa'};
            for ii = 1:3
                fidOrder(ii) = find(ismember(lower(shape.fid.label),targetOrder{ii}));
            end
            fid.fid.label = targetOrder;
            fid.fid.pnt = shape.fid.pos(fidOrder,:);
            
        else % coordsys.json fallback
            
            try
                coord = spm_load(S.coordsystem);
            catch
                coord = spm_load(coordFile);
            end
            
            % These HeadCoilCoordinates (nas,lpa,rpa) are same space as
            % sensors positions and specified in the coordsystem.json file
            fiMat(1,:) = coord.HeadCoilCoordinates.coil1;
            fiMat(2,:) = coord.HeadCoilCoordinates.coil2;
            fiMat(3,:) = coord.HeadCoilCoordinates.coil3;
            
            fid.fid.label = {'nas', 'lpa', 'rpa'}';
            fid.fid.pnt = fiMat;
            
            fid.pnt = []; % headshape field that is left blank,
            
        end
        
        % Use SPM Template brain template
        M.fid.label = {'nas', 'lpa', 'rpa'}';
        M.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:);
        M.fid.pos = []; % headshape field that is left blank (GRB)
        M.pnt = D.inv{1}.mesh.fid.pnt;
    catch
        subjectNoStruct = 0;
        warning(['Could not load coordinate system - assuming sensor ',...
            'positions are in the same coordinate as SPM canonical brain']);
    end
end

if(template && ~subjectNoStruct) %make
    fid.fid.label = {'nas', 'lpa', 'rpa'}';
    fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:);
    fid.pos= []; % headshape field that is left blank (GRB)
    M = fid;
    M.pnt = D.inv{1}.mesh.fid.pnt;
end

%- Coregistration
%--------------------------------------------------------------------------
if(forward)
    if(subjectSource||template)
        D = fiducials(D, fid);
        save(D);
        f=fiducials(D);
        if ~isempty(f.pnt)
            D = spm_eeg_inv_datareg_ui(D,1,f,M,1);
        else
            f.pnt = zeros(0,3);
            D = spm_eeg_inv_datareg_ui(D,1,f,M,0);
        end
    end
end
%- Foward  model specification
%--------------------------------------------------------------------------
if forward
    D.inv{1}.forward.voltype = S.voltype;
    D = spm_eeg_inv_forward(D);
    spm_eeg_inv_checkforward(D,1,1);
end
save(D);

%- Create lead fields
%--------------------------------------------------------------------------
if(forward)
    D.inv{1}.forward.voltype = S.voltype;
    D = spm_eeg_inv_forward(D);
    nverts = length(D.inv{1}.forward.mesh.vert);
    if(S.lead)
        [L,D] = spm_eeg_lgainmat(D,1:nverts);
    end
    spm_eeg_inv_checkforward(D,1,1);
    save(D);
end
fprintf('%-40s: %30s\n','Completed',spm('time'));


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
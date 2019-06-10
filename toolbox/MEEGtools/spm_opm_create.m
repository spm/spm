function [D,L] = spm_opm_create(S)
% Read or simulate magnetometer data and optionally set up forward model
% FORMAT D = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
% SENSOR LEVEL INFO
%   S.data          - filepath to .bin file        - Default: Simulates data 
%   S.channels      - channels.tsv file            - Default: REQUIRED  
%   S.fs            - Sampling frequency (Hz)      - Default: REQUIRED if S.meg is empty
%   S.meg           - meg.json file                - Default: REQUIRED if S.fs is empty
% SIMULATION
%   S.wholehead     - whole head coverage flag     - Deafult: 0
%   S.space         - space between sensors(mm)    - Default: 25
%   S.offset        - scalp to sensor distance(mm) - Default: 6.5
%   S.nSamples      - number of samples            - Default: 1000
% SOURCE LEVEL INFO
%   S.coordsystem   - coordsystem.json file        - Default: 
%   S.positions     - positions.tsv file           - Default:
%   S.sMRI          - Filepath to  MRI file        - Default: uses template
%   S.cortex        - Custom cortical mesh         - Default: Use inverse normalised cortical mesh
%   S.scalp         - Custom scalp mesh            - Default: Use inverse normalised scalp mesh
%   S.oskull        - Custom outer skull mesh      - Default: Use inverse normalised outer skull mesh
%   S.iskull        - Custom inner skull mesh      - Default: Use inverse normalised inner skull mesh
%   S.voltype       - Volume conducter Model type  - Default: 'Single Shell'
%   S.meshres       - mesh resolution(1,2,3)       - Default: 1
%   S.lead          - flag to compute lead field   - Default: 0
% Output:
%  D           - MEEG object (also written to disk)
%  L           - Lead field (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_create.m 7612 2019-06-10 16:24:05Z tim $
spm('FnBanner', mfilename);

%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'voltype'),     S.voltype = 'Single Shell'; end
if ~isfield(S, 'meshres'),     S.meshres = 1; end
if ~isfield(S, 'scalp'),       S.scalp = []; end
if ~isfield(S, 'cortex'),      S.cortex = []; end
if ~isfield(S, 'iskull'),      S.iskull = []; end
if ~isfield(S, 'oskull'),      S.oskull = []; end
if ~isfield(S, 'lead'),        S.lead = 0; end
if ~isfield(S, 'fs'),          S.fs   = 1000; end
if ~isfield(S, 'nSamples'),    S.nSamples   = 1000; end
if ~isfield(S, 'space'),       S.space  = 25; end
if ~isfield(S, 'offset'),      S.offset  = 6.5; end
if ~isfield(S, 'data'),        S.data = zeros(1,S.nSamples); end
if ~isfield(S, 'wholehead'),   S.wholehead = 1; end
if ~isfield(S, 'fname'),       S.fname = 'sim_opm'; end
     
%-Simulate or read
%--------------------------------------------------------------------------
readData= isa(S.data,'char');

if readData
    [a,b,~] = fileparts(S.data);
    
    %- Check for channel Info
    %----------------------------------------------------------------------
    if ~isfield(S,'channels')
        base = strsplit(b,'meg');
        chanFile= fullfile(a,[base{1},'channels.tsv']);
        if exist(chanFile,'file')==2
            channels = spm_load(chanFile);
        else
            error('A channels.tsv file is necessray to read the data');
        end
    else
        if exist(S.channels,'file')==2
            channels = spm_load(S.channels);
        else
            error('The file in S.channels deos not exist');
        end
    end
    
    %- Check for MEG Info
    %----------------------------------------------------------------------
    if ~isfield(S,'meg')
        base = strsplit(b,'meg');
        megFile= fullfile(a,[base{1},'meg.json']);
        if exist(megFile,'file')==2
            meg = spm_load(megFile);
        else
            if isfield(S,'fs')
                meg=[];
                meg.SamplingFrequecny=S.fs;
            else
                error('Either S.fs or S.meg should be specified')
            end
        end
    else
        if exist(S.meg,'file')==2
            meg = spm_load(S.meg);
        else
            error('The file in S.meg deos not exist');
        end
    end
    
    %- Read Binary File
    %----------------------------------------------------------------------
    if exist(S.data,'file')==2
        dat = fopen(S.data);
        bytes = fread(dat,Inf,'double',0,'b');
        fclose(dat);
        nc = size(channels.name,1);
        S.data = reshape(bytes,nc,length(bytes)/nc);
    else
        error('The file in S.data deos not exist');
    end
    
    %- Forward model Check
    %----------------------------------------------------------------------
    if ~isfield(S, 'sMRI')
        forward=0;
    elseif ~exist(S.sMRI,'file')==2
        error('S.sMRI file does not exist.')
    else 
        forward =1;
        template=0;
    end
else
    %- Simulation Defaults
    %----------------------------------------------------------------------
    
    a=pwd();
    b=S.fname;
    
    args=[];
    args.base='Chan';
    args.n= size(S.data,1);
    labs= spm_create_labels(args);
    channels = [];
    channels.name=labs;
    channels.type=repmat({'MEG'},size(S.data,1),1);
    channels.units= repmat({'fT'},size(S.data,1),1);
    meg=[];
    meg.SamplingFrequency=S.fs;
    
    forward=1;
    if ~isfield(S, 'sMRI')
        S.sMRI   = 1;
        template = 1;
    else
        template = 0;
    end
end 

%- Create SPM object of simulated or real data
%--------------------------------------------------------------------------
Dtemp = meeg(size(S.data,1),size(S.data,2),size(S.data,3));
Dtemp = fsample(Dtemp,meg.SamplingFrequency);
Dtemp = fname(Dtemp,[b,'.mat']);
Dtemp = path(Dtemp,a);
Dtemp = chanlabels(Dtemp,1:size(Dtemp,1),channels.name);
Dtemp = units(Dtemp,1:size(Dtemp,1),channels.units);
Dtemp = chantype(Dtemp,1:size(Dtemp,1),channels.type);

%- Overwrite and Save
%--------------------------------------------------------------------------
ma = fullfile(a,[b,'.mat']);
da = fullfile(a,[b,'.dat']);

ae = exist(fname(Dtemp),'file')==2;
if(ae)
    delete(ma);
    delete(da);
end
Dtemp.save();
% create data file and insert data
D= blank(Dtemp,[b,'.dat']);
dim=size(D);
D(1:dim(1),1:dim(2),1:dim(3)) = S.data;
D.save();

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

%- Create the Sensor Array
%--------------------------------------------------------------------------
if forward
    if ~readData
        if isfield(S, 'positions')
            posOri = spm_load(S.positions,'',true);
            pos = [posOri.Px,posOri.Py,posOri.Pz];
            ori = [posOri.Ox,posOri.Oy,posOri.Oz];
            
        else
            % if no postions and orientations provided then create them
            args = [];
            args.D =D;
            args.offset = S.offset;
            args.space = S.space;
            args.wholehead = S.wholehead;
            [pos,ori] = opm_createSensorArray(args);
        end
        
        nSensors=size(pos,1);
        args=[];
        args.base='Chan';
        args.n= nSensors;
        labs= spm_create_labels(args);
        channels = [];
        channels.name=labs;
        channels.type=repmat({'MEG'},nSensors,1);
        channels.units= repmat({'fT'},nSensors,1);
        D = clone(D,fnamedat(D),[nSensors,S.nSamples,1],1);
        D = chanlabels(D,1:size(D,1),channels.name);
        D = units(D,1:size(D,1),channels.units);
        D = chantype(D,1:size(D,1),channels.type);
        cl=chanlabels(D)';
    else
        if ~isfield(S,'positions')
            base = strsplit(b,'meg');
            posFile= spm_select('FPList',a,[base{1},'positions.tsv']);
            posOri = spm_load(posFile);
            pos = [posOri.Px,posOri.Py,posOri.Pz];
            ori = [posOri.Ox,posOri.Oy,posOri.Oz];
            cl=posOri.name;
        else
            if exist(S.positions,'file')==2
               posOri = spm_load(S.positions);
                pos = [posOri.Px,posOri.Py,posOri.Pz];
                ori = [posOri.Ox,posOri.Oy,posOri.Oz];
                cl=posOri.name;
            else
                error('The file in S.meg deos not exist');
            end
        end
    end
end

%-Place Sensors  in object
%--------------------------------------------------------------------------
if forward
    grad= [];
    grad.label = cl;
    grad.coilpos = pos;
    grad.coilori = ori;
    grad.tra = eye(numel(grad.label));
    grad.chanunit = repmat({'T'}, numel(grad.label), 1);
    grad.chantype= 'MEG';
    grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'mm');
    D = sensors(D, 'MEG', grad);
    save(D);
end


%- fiducial settings
%--------------------------------------------------------------------------
if forward
    if readData
        if ~isfield(S,'coordsystem')
            base = strsplit(b,'meg');
            coordFile= fullfile(a,[base{1},'coordsystem.json']);
            if exist(coordFile,'file')==2
                coord = spm_load(coordFile);
            else
                error('A coordsystem.json file is required for forward modelling')
            end
        else
            if exist(S.coordsystem,'file')==2
                coord = spm_load(S.coordsystem);
            else
                error('The file in S.coordsystem deos not exist');
            end
        end
        
        fid.fid.label = fieldnames(coord.HeadCoilCoordinates);
        f1=coord.HeadCoilCoordinates.coil1;
        f2=coord.HeadCoilCoordinates.coil2;
        f3=coord.HeadCoilCoordinates.coil3;
        fid.fid.pnt =[f1';f2';f3'];
        fid.pos= []; % headshape field that is left blank (GRB)
       
        m1=coord.AnatomicalLandmarkCoordinates.coil1;
        m2=coord.AnatomicalLandmarkCoordinates.coil2;
        m3=coord.AnatomicalLandmarkCoordinates.coil3;
        M = fid;
        M.fid.pnt=[m1';m2';m3'];
        M.pnt = D.inv{1}.mesh.fid.pnt;
        
    elseif(template)
        fid.fid.label = {'nas', 'lpa', 'rpa'}';
        fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:);
        fid.pos= []; % headshape field that is left blank (GRB)
        M = fid;
        M.pnt = D.inv{1}.mesh.fid.pnt;
    else
        fid.fid.label = {'nas', 'lpa', 'rpa'}';
        fid.fid.pnt = [0 0 0; -1 0 0; 1 0 0];
        fid.pos= [];
        M = fid;
        M.pnt = D.inv{1}.mesh.fid.pnt;
    end
    
    D = fiducials(D, fid);
    save(D);
    f=fiducials(D);
    f.pnt =zeros(0,3);
    D = spm_eeg_inv_datareg_ui(D,1,f,M,0);
end

%- 2D view based on mean orientation of sensors 
%--------------------------------------------------------------------------
if(forward)
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
%- Foward  model specification
%--------------------------------------------------------------------------
if forward
    D.inv{1}.forward.voltype = S.voltype;
    D = spm_eeg_inv_forward(D);
    nverts = length(D.inv{1}.forward.mesh.vert);
    if(S.lead)
    [L,D] = spm_eeg_lgainmat(D,1:nverts);
    end
    spm_eeg_inv_checkforward(D,1,1);
end
save(D);
fprintf('%-40s: %30s\n','Completed',spm('time'));    


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Create Sensor Array                                 %                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos,ori] = opm_createSensorArray(S)
% Given a scalp surface even samples the surface with sensors 
% 
% FORMAT [pos,ori] = spm_opm_createSensorArray(S)
%   S               - input structure
% Fields of S:
%   S.D            - SPM M/EEG object with surface meshes 
%   S.offset       - distance to place sensors(mm) from scalp surface
%   S.space        - distance between sensors(mm) 
%   S.wholehead    - boolean: Should whole scalp surface should be covered?
% _________________________________________________________________________

% Args
%--------------------------------------------------------------------------
D = S.D;
wholehead = S.wholehead;

% Meshes
%--------------------------------------------------------------------------
scalp = gifti(D.inv{1}.mesh.tess_scalp);
cortex = gifti(D.inv{1}.mesh.tess_ctx);
lp = min(cortex.vertices(:,3));

% create convex hull of mesh
%--------------------------------------------------------------------------
[~,V1] = convhull(double(scalp.vertices));

% outward facing vertx normals of convex hull
%--------------------------------------------------------------------------
[Nv,~] = spm_mesh_normals(scalp,true);
ns = Nv;
vs = scalp.vertices;
cog1=mean(vs);

for i = 1:length(ns)
    ad = vs(i,:) + 5*ns(i,:);
    subtrac = vs(i,:) - 5*ns(i,:);
    
    d1 = sum((cog1 - ad).^2);
    d2 = sum((cog1 - subtrac).^2);
    
    if(d2>d1)
        ns(i,:) = -ns(i,:);
    else
        ns(i,:) =ns(i,:);
    end
end

% add an offset to the convex hull
%--------------------------------------------------------------------------
offVertices=vs+ns*S.offset;
[~,V2] = convhull(double(offVertices));


% Expand solution space by ratio of Volumes
%--------------------------------------------------------------------------
T= eye(4)*(V2/V1)^(1/3);
T(4,4)=1;
scalp = spm_mesh_transform(scalp,T);

% Translate solution space so centre of gravities overlap
%--------------------------------------------------------------------------

cog2 = mean(scalp.vertices);
T= eye(4);
T(1:3,4)=cog1-cog2+0;
scalp = spm_mesh_transform(scalp,T);

% Create the sensor array 
%--------------------------------------------------------------------------
args= [];
args.space=S.space;
args.g=scalp;
[pos, ~, ~] = spm_mesh_pack_points(args);


% get orientation of scalp
%--------------------------------------------------------------------------
[~,Nf] = spm_mesh_normals(scalp,true);
cog=mean(scalp.vertices);

v=scalp.vertices;
f=scalp.faces;
faceP=zeros(size(f,1),3);
for i =1:length(faceP)
    whichVerts= f(i,:);
    verts= v(whichVerts,:);
    faceP(i,:)=mean(verts);
end

% make all normals point outwards
%--------------------------------------------------------------------------
for i = 1:length(Nf)
    ad = faceP(i,:) + 5*Nf(i,:);
    subtrac = faceP(i,:) - 5*Nf(i,:);
    
    d1 = sum((cog - ad).^2);
    d2 = sum((cog - subtrac).^2);
    
    if(d2>d1)
        Nf(i,:) = -Nf(i,:);
    else
        Nf(i,:) =Nf(i,:);
    end
end

% assign orientatation to sensors of closest face normal
%--------------------------------------------------------------------------
  ori=zeros(size(pos));
for i =1:length(ori)
    tmp=pos(i,:);
    di=sqrt(sum(bsxfun(@minus,faceP,tmp).^2,2));
    [~,indmin]=min(di);
    ori(i,:)=-Nf(indmin,:);
end

% Check if wholehead is requested
%--------------------------------------------------------------------------
if(~wholehead)
    C= pos(:,3)>(lp);
    pos= pos(C,:);
    ori = ori(C,:);
end

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
 
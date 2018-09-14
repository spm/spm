function [D,L] = spm_opm_create(S)
% Create a valid M/EEG object from raw data or simulate magnetometer data
% FORMAT D = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
%   S.data          - 2/3 dimensaional array       - Default: empty dataset 
%   S.pinout        - filepath to pinout file      - Default: no labels
%   S.fname         - filename for  dataset        - Default: 'OPM.dat'
%   S.fs            - Sampling rate of dataset     - Default: 1000
%   S.scale         - scale factor (to fT)         - Default: 1
%   S.trig          - trigger matrix               - Default: no triggers 
%   S.pinout        - Mapping of labels to S.data  - Default: generate labels 
%   S.sensorsUsed   - organisation of S.data       - Default: all of S.data is considered MEG data 
%   S.sMRI          - Filepath to  MRI file        - Default: uses template
%   S.cortex        - Custom cortical mesh         - Default: Use inverse normalised cortical mesh
%   S.scalp         - Custom scalp mesh            - Default: Use inverse normalised scalp mesh
%   S.oskull        - Custom outer skull mesh      - Default: Use inverse normalised outer skull mesh
%   S.iskull        - Custom inner skull mesh      - Default: Use inverse normalised inner skull mesh
%   S.voltype       - Volume conducter Model type  - Default: 'Single Shell'
%   S.meshres       - mesh resolution(1,2,3)       - Default: 1
%   S.fid           - 3 x 3 Fiducial  matrix       - Default: [0 0 0; -1 0 0; 1 0 0]
%   S.wholehead     - whole head coverage flag     - Deafult: 0
%   S.space         - space between sensors(mm)    - Default: 25
%   S.offset        - scalp to sensor distance(mm) - Default: 6.5
%   S.nSamples      - number of samples            - Default: 1000
%   S.lead          - flag to compute lead field   - Default: 0
% Output:
%  D           - MEEG object (also written to disk)
%  L           - Lead field (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_create.m 7418 2018-09-14 13:16:33Z tim $


%-Initialise
%--------------------------------------------------------------------------

% sensor level data with(out) forward model
if (~isfield(S, 'sMRI') && isfield(S, 'data'))
    if isfield(S, 'pos')
        error('if S.pos is not empty a structural  MRI is required')
    end
    forward = 0;
else
    forward = 1;
end

% template based simulation
if ~isfield(S, 'sMRI') 
    S.sMRI   = 1;
    template = 1;
else
    template = 0;
end

% labeled data
if (isfield(S,'pinout') && isfield(S,'data'))    
    % read the pinout and opm2cast file
    pinout = spm_load(S.pinout,'',false);
    sensorsUsed = spm_load(S.sensorsUsed,'',false);  
    labeledData = 1;
else
    labeledData = 0;
end


%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'fname'),       S.fname   = 'OPM.dat'; end
if ~isfield(S, 'fs'),          S.fs   = 1000; end
if ~isfield(S, 'nSamples'),    S.nSamples   = 1000; end
if ~isfield(S, 'space'),       S.space  = 25; end
if ~isfield(S, 'offset'),      S.offset  = 6.5; end
if ~isfield(S, 'fid'),         S.fid  = [0 0 0; -1 0 0; 1 0 0]; end
if ~isfield(S, 'voltype'),     S.voltype = 'Single Shell'; end
if ~isfield(S, 'meshres'),     S.meshres = 1; end
if ~isfield(S, 'wholehead'),   S.wholehead = 0; end
if ~isfield(S, 'scale'),       S.scale = 1; end
if ~isfield(S, 'data'),        S.data = zeros(1,S.nSamples); end
if ~isfield(S, 'scalp'),       S.scalp = []; end
if ~isfield(S, 'cortex'),      S.cortex = []; end
if ~isfield(S, 'iskull'),      S.iskull = []; end
if ~isfield(S, 'oskull'),      S.oskull = []; end
if ~isfield(S, 'lead'),        S.lead = 0; end

%-File Management
%--------------------------------------------------------------------------
[a,b]= fileparts(S.fname);
ma = fullfile(a,[b,'.mat']);
da = fullfile(a,[b,'.dat']);

if exist(ma,'file')==2
    delete(ma);
end
if exist(da,'file')==2
    delete(da);
end

%- Select Channels
%--------------------------------------------------------------------------
% labelled data of different types(MEG, REF)
if labeledData
    used = zeros(size(sensorsUsed.Var1,1),1);
    for i = 1:length(used)
        used(i) = find(strcmp(sensorsUsed.Var1{i},pinout.Var2));
    end
    matPos= pinout.Var1(used);
    labs = pinout.Var2(used);
    refInd = find(strcmp('REF',sensorsUsed.Var2));
    nMeg = length(used)-length(refInd);
    S.data = S.data(matPos,:);
    megInd= setdiff(1:size(used,1),refInd)';
else
    refInd=[];
    nMeg = size(S.data,1);
    megInd = 1:nMeg;
end
%- Account for triggers
%--------------------------------------------------------------------------
if isfield(S ,'trig')
    binTrig = zeros(size(S.trig));
    S.data = [S.data;binTrig];
    for j =1:size(S.trig,1)
        binTrig(j,:)= S.trig(j,:);
    end
    St = [];
    St.base='TRIG';
    St.n=size(binTrig,1);
    triglabs = spm_create_labels(St);
else
    binTrig = [];
    triglabs={};
end
%- Account for PHYS type
%--------------------------------------------------------------------------
if isfield(S ,'other')
    other = zeros(size(S.other));
    S.data = [S.data;other];
    St = [];
    St.base='PHYS';
    St.n=size(S.other,1);
    physlabs = spm_create_labels(St);
else
    physlabs={};
end
%- Account for different channel types
%--------------------------------------------------------------------------
% initially say its all MEG data
cType = {};
cType{1} = 'MEG';
cType = repmat(cType,size(S.data,1),1);

% then add REF labels if they exist
for i = 1:length(refInd)
    cType{refInd(i),:} = 'REF';
end
% then add TRIG type if they exist
for i = 0:(size(binTrig,1)-1)
    cType{(end-i)} ='TRIG';
    trigInd = find(strcmp('TRIG',cType));
end

% then add PHYS type if they exist
if(isfield(S,'other'))
stPhys=size(cType,1)-size(binTrig,1)-size(S.other,1)+1;
endPhys=size(cType,1)-size(binTrig,1);
for i = stPhys:endPhys
    cType{i,:} ='PHYS';
end
physInd = find(strcmp('PHYS',cType));
end
%-Sensor Level Data
%--------------------------------------------------------------------------
D = opm_convert(S.data,S.fname,S.fs,S.scale);

% initially D is filled with zeros in trigger channel. This prevents
% scaling of the triggers
if isfield(S ,'trig')
    D(trigInd,:,:)= binTrig;
    save(D);
end

if isfield(S ,'other')
    D(physInd,:,:)= S.other;
    save(D);
end
%-Scalp Extraction
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
    % if user supplies postions and orientations
    if isfield(S, 'pos')
        posOri = spm_load(S.pos,'',false);
        
        % if its labelled data subset the postions and orientaitons
        if(labeledData)
            megInd = find(strcmp('MEG',cType));
            slots = sensorsUsed.Var2(megInd)'; 
            posSlots = posOri.Var7;
            slotind = zeros(length(slots),1);
            for j = 1:length(slots)
                slotind(j)= find(strcmp(slots{j},posSlots));
            end
            px = posOri.Var1(slotind);
            py = posOri.Var2(slotind);
            pz = posOri.Var3(slotind);
            pos = [px,py,pz];
            ox = posOri.Var4(slotind);
            oy = posOri.Var5(slotind);
            oz = posOri.Var6(slotind);
            ori = [ox,oy,oz];
            nSensors = size(pos,1);
            
            % if not labeled data then don't subset
        else
            pos = [posOri.Var1,posOri.Var2,posOri.Var3];
            ori = [posOri.Var4,posOri.Var5,posOri.Var6];
            nSensors = size(pos,1);
            D = clone(D,S.fname,[nSensors,S.nSamples,1],1);
            megInd = 1:nSensors;
        end
        
        % if no postions and orientations provided then create them
    else
        args = [];
        args.D =D;
        args.offset = S.offset;
        args.space = S.space;
        args.wholehead = S.wholehead;
        [pos,ori] = opm_createSensorArray(args);
        nSensors = length(pos);
        D = clone(D,S.fname,[nSensors,S.nSamples,1],1);
        megInd = 1:nSensors;
    end
end

%-units, labels, and types
%--------------------------------------------------------------------------
nSensors = size(D,1);

% update channel type (sensor number changes when simulating)
if length(cType)<nSensors
    cType = repmat(cType,nSensors,1);
end

% if we don't have labelled data then create labels
if ~labeledData
    St= [];
    St.base ='Sensor';
    St.n = nSensors;
    labs = spm_create_labels(St);
end

% now we know how many sensors we have we can set units labels and types
labs = [labs;physlabs';triglabs'];
D = units(D,[megInd;refInd],'fT');
D = chantype(D,1:nSensors,cType);
D = chanlabels(D,1:size(D,1),labs');
save(D);

%-Place Sensors  in object
%--------------------------------------------------------------------------
if forward
    grad= [];
    grad.label = {labs{megInd}}';
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
% because we're simulating we already know where things are so no need to
% supply realistic values unless required for compatibility with other
% code. If using a template fiducials are adjusted to be the MNI
% fiducials. MUST UPDATE THIS TO PUT IN VAGUELY SENSIBLE FIDUCIALS BY
% TRANSFORMING MNI FIDUCIALS INTO THE INDIVIDUAL SPACE...

if forward
    fid.fid.label = {'nas', 'lpa', 'rpa'}';
    if(template)
        fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:);
    else
        fid.fid.pnt = S.fid;
    end
    fid.pos= []; % headshape field that is left blank (GRB)
    D = fiducials(D, fid);
    save(D);
end
%- Coregister
%--------------------------------------------------------------------------
% We already have the meshes and sensors in same space but we need to add
% the datareg field so SPM won't get confused at the forward model stage
if forward
    f = fiducials(D);
    f.pnt =zeros(0,3);
    M = f;
    M.pnt = D.inv{1}.mesh.fid.pnt;
    D = spm_eeg_inv_datareg_ui(D,1,f,M,0);
end

%- 2D view based on mean orientation of sensors 
%--------------------------------------------------------------------------
if(forward)
    n1=mean(grad.coilori); n1= n1./sqrt(dot(n1,n1));
    t1=cross(n1,[0 0 1]);
    t2=cross(t1,n1);
    
    for i=1:size(grad.coilpos,1)
        pos2d(i,1)=dot(grad.coilpos(i,:),t1);
        pos2d(i,2)=dot(grad.coilpos(i,:),t2);
    end
    
    
    args=[];
    args.D=D;
    args.xy= pos2d';
    args.label=grad.label;
    args.task='setcoor2d';
    D=spm_eeg_prep(args);
    D.save;
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
    spm_eeg_inv_checkforward(D,1,1)
end
save(D)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           OPM CONVERT                                   %                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dout = opm_convert(array,fnamedat,fs,scale)
% Convert array into SPM MEEG object
% FORMAT Dout = spm_opm_convert(array,fnamedat,fs,scale)
%
% array       - numeric Array of 2,3,4 dimensions(channels,time,trials)
% fnamedat    - string specifing output path of object(include extension .dat)
% fs          - sampling frequency
% scale       - scale factor to convert to fT [default: 1]
%__________________________________________________________________________
% Copyright (C) 2017 Tim Tierney

% Tim Tierney

% determine output filename
[a, b] = fileparts(fnamedat);
outMat = fullfile(a,[b,'.mat']);

% if scale is not supplied set a default  of 1 
if nargin < 4
    scale = 1;
end
    
array = array.*scale;

% find number of dimensions and decide what to do based on result
dim = size(array);
L   = length(dim);

if L > 4
    % throw error for having unsupported number of dimensions
    error('Array should not have more than 4 dimensions.');
elseif L == 4
    % create MEG object with appropriate size
    Dout = meeg(dim(1),dim(2),dim(3),dim(4));
    
    % make it a blank object
    Dout= blank(Dout,fnamedat);
    
    % set the smaple rate 
    Dout = Dout.fsample(fs);
    
    % fill in data with supplied array
    Dout(1:dim(1),1:dim(2),1:dim(3),1:dim(4)) = array;
elseif L == 3
    % same with less dimensions
    Dout = meeg(dim(1),dim(2),dim(3));
    Dout= blank(Dout,fnamedat);
    Dout = Dout.fsample(fs);
    Dout(1:dim(1),1:dim(2),1:dim(3)) = array;
    
elseif L == 2 
    %same with less dimensions
    Dout = meeg(dim(1),dim(2),1);
    Dout= blank(Dout,fnamedat);
    Dout = Dout.fsample(fs);
    Dout(1:dim(1),1:dim(2),1) = array;
else
    % throw exception if I really don't know what to do
    error('Array must have between 2 and  4 dimensions.');
end

% set the filename and save
Dout = fname(Dout,outMat);
Dout = path(Dout,a);
Dout.save;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Grid Surface Intersection                           %                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outs = poiGridSurface(surface,space,lowerThresh,upperThresh)
% Find where a grid intersects a surface (ideally convex)

% Create Grid 
%--------------------------------------------------------------------------
bbMin = min(surface.vertices);
bbMax = max(surface.vertices);
x = (bbMin(1)-space):space:(bbMax(1)+space);
y = (bbMin(2)-space):space:(bbMax(2)+space);
z = (bbMin(3)-space):space:(bbMax(3)+space);
[X,Y,Z] = ndgrid(x,y,z);
grid   = [X(:) Y(:) Z(:)];

% Keep only outer layer of grid
%--------------------------------------------------------------------------
mi = min(grid);
ma = max(grid);

miX = grid(find(abs(grid(:,1) - mi(1))<=.01),:);
miY = grid(find(abs(grid(:,2) - mi(2))<=.01),:);
miZ = grid(find(abs(grid(:,3) - mi(3))<=.01),:);

maX = grid(find(abs(grid(:,1) - ma(1))<=.01),:);
maY = grid(find(abs(grid(:,2) - ma(2))<=.01),:);
maZ = grid(find(abs(grid(:,3) - ma(3))<=.01),:);




% ny
%--------------------------------------------------------------------------
outs = zeros(length(maY),3);

for j = 1:length(maY)
Ia = maY(j,:);
Ib = maY(j,:)+200.*[0,-1,0];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
m(:,2) = verts(2,:)-verts(1,:);
m(:,3) = verts(3,:)-verts(1,:);
Y = (Ia-verts(1,:))';
tuv(i,:) = m\Y;
end

a = tuv(:,1)< 1 & tuv(:,1)>0;
b = tuv(:,2)< 1 & tuv(:,2)>0;
c = tuv(:,3)< 1 & tuv(:,3)>0;
d = (tuv(:,2)+tuv(:,3))<=1;
e = a&b&c&d;
mul =min(tuv(e,1));

if(~isempty(mul))
outs(j,:) = Ia+(Ib-Ia)*mul;
end

end
ind = sum(abs(outs),2)>0;
outsny = outs(ind,:);
disp('Negative Y')

% py
%--------------------------------------------------------------------------
outs = zeros(length(miY),3);

for j = 1:length(miY)
Ia = miY(j,:);
Ib = miY(j,:)+200.*[0,1,0];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
m(:,2) = verts(2,:)-verts(1,:);
m(:,3) = verts(3,:)-verts(1,:);
Y = (Ia-verts(1,:))';
tuv(i,:) = m\Y;
end

a = tuv(:,1)< 1 & tuv(:,1)>0;
b = tuv(:,2)< 1 & tuv(:,2)>0;
c = tuv(:,3)< 1 & tuv(:,3)>0;
d = (tuv(:,2)+tuv(:,3))<=1;
e = a&b&c&d;
mul =min(tuv(e,1));

if(~isempty(mul))
outs(j,:) = Ia+(Ib-Ia)*mul;
end

end
ind = sum(abs(outs),2)>0;
outspy = outs(ind,:);
disp('Positive Y')

% nx
%--------------------------------------------------------------------------
outs = zeros(length(maX),3);

for j = 1:length(maX)
Ia = maX(j,:);
Ib = maX(j,:)+200.*[-1,0,0];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
m(:,2) = verts(2,:)-verts(1,:);
m(:,3) = verts(3,:)-verts(1,:);
Y = (Ia-verts(1,:))';
tuv(i,:) = m\Y;
end

a = tuv(:,1)< 1 & tuv(:,1)>0;
b = tuv(:,2)< 1 & tuv(:,2)>0;
c = tuv(:,3)< 1 & tuv(:,3)>0;
d = (tuv(:,2)+tuv(:,3))<=1;
e = a&b&c&d;
mul =min(tuv(e,1));

if(~isempty(mul))
outs(j,:) = Ia+(Ib-Ia)*mul;
end

end
ind = sum(abs(outs),2)>0;
outsnx = outs(ind,:);
disp('Negative X');

% px
%--------------------------------------------------------------------------
outs = zeros(length(miX),3);

for j = 1:length(miX)
Ia = miX(j,:);
Ib = miX(j,:)+200.*[1,0,0];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
m(:,2) = verts(2,:)-verts(1,:);
m(:,3) = verts(3,:)-verts(1,:);
Y = (Ia-verts(1,:))';
tuv(i,:) = m\Y;
end

a = tuv(:,1)< 1 & tuv(:,1)>0;
b = tuv(:,2)< 1 & tuv(:,2)>0;
c = tuv(:,3)< 1 & tuv(:,3)>0;
d = (tuv(:,2)+tuv(:,3))<=1;
e = a&b&c&d;
mul =min(tuv(e,1));

if(~isempty(mul))
outs(j,:) = Ia+(Ib-Ia)*mul;
end

end
ind = sum(abs(outs),2)>0;
outspx = outs(ind,:);
disp('Postive X')

% nz
%--------------------------------------------------------------------------
outs = zeros(length(maZ),3);

for j = 1:length(maZ)
Ia = maZ(j,:);
Ib = maZ(j,:)+200.*[0,0,-1];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
m(:,2) = verts(2,:)-verts(1,:);
m(:,3) = verts(3,:)-verts(1,:);
Y = (Ia-verts(1,:))';
tuv(i,:) = m\Y;
end

a = tuv(:,1)< 1 & tuv(:,1)>0;
b = tuv(:,2)< 1 & tuv(:,2)>0;
c = tuv(:,3)< 1 & tuv(:,3)>0;
d = (tuv(:,2)+tuv(:,3))<=1;
e = a&b&c&d;
mul =min(tuv(e,1));

if(~isempty(mul))
outs(j,:) = Ia+(Ib-Ia)*mul;
end

end
ind = sum(abs(outs),2)>0;
outsnz = outs(ind,:);

disp('Negative Z');

% pz
%--------------------------------------------------------------------------
outs = zeros(length(miZ),3);

for j = 1:length(miZ)
Ia = miZ(j,:);
Ib = miZ(j,:)+200.*[0,0,1];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
m(:,2) = verts(2,:)-verts(1,:);
m(:,3) = verts(3,:)-verts(1,:);
Y = (Ia-verts(1,:))';
tuv(i,:) = m\Y;
end

a = tuv(:,1)< 1 & tuv(:,1)>0;
b = tuv(:,2)< 1 & tuv(:,2)>0;
c = tuv(:,3)< 1 & tuv(:,3)>0;
d = (tuv(:,2)+tuv(:,3))<=1;
e = a&b&c&d;
mul =min(tuv(e,1));

if(~isempty(mul))
outs(j,:) = Ia+(Ib-Ia)*mul;
end

end
ind = sum(abs(outs),2)>0;
outspz = outs(ind,:);

disp('Positive Z');

% Initial Points
%--------------------------------------------------------------------------
uspz = unique(round(outsnz,3),'rows');
% nn= [];
% rem = zeros(length(uspz),1);
% dist = squareform(pdist(uspz));
% 
% for i =1:length(uspz)
%     sdist= sort(dist(:,i));
%     nn= min(sdist(2:4));
%     if(nn<upperThresh)
%         rem(i)=1;
%     else
%         rem(i)=0;
%     end
%         
% end
% uspz = uspz(boolean(rem),:);

% Constrained adding of points
%--------------------------------------------------------------------------
other = [unique(round(outsnx,3),'rows');...
         unique(round(outspy,3),'rows');...
         unique(round(outsny,3),'rows');...
         unique(round(outspz,3),'rows');...
         unique(round(outspx,3),'rows');...
         ];
    us = uspz;
    
for i = 1:length(other)
ot = other(i,:);
di = sqrt(sum((bsxfun(@minus,us,ot)).^2,2));
if(all(di>lowerThresh))
   us = [us;ot];
   
% rem = zeros(length(us),1);
% dist = squareform(pdist(us));
% for i =1:length(us)
%     sdist= sort(dist(:,i));
%     nn= min(sdist(2:6));
%     if(nn>upperThresh)
%         rem(i)=0;
%     else
%         rem(i)=1;
%     end
%         
% end 
end

end

% return
%--------------------------------------------------------------------------     
outs = us;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Line Surface Intersection                           %                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outs = poiLineSurface(surface,points,normals,dist)
% Find where a point intersects a surface 

%Point of intersection
%--------------------------------------------------------------------------
outs = zeros(length(points),3);

% compute all face centres
posi = zeros(size(surface.faces));
for i = 1:length(posi)
    whichVerts = surface.faces(i,:);
    cos = surface.vertices(whichVerts,:);
    posi(i,:) = mean(cos);
end


for j = 1:length(points)
    Ia = points(j,:);
    Ib = points(j,:)+dist.*normals(j,:);
    
    m = zeros(3,3);
    m(:,1) = Ia-Ib;
    
    di = sqrt(sum(bsxfun(@minus,posi,Ia).^2,2));
    faceInd = find(di<dist);
    
    
    smallFaces = surface.faces(faceInd,:);
    nfaces = length(smallFaces);
    
    tuv = zeros(nfaces,3);
    verts = double(zeros(3,3));
    Y = zeros(3,1);
    
    for i = 1:nfaces
        verts = surface.vertices(smallFaces(i,:),:);
        m(:,2) = verts(2,:)-verts(1,:);
        m(:,3) = verts(3,:)-verts(1,:);
        Y = (Ia-verts(1,:))';
        tuv(i,:) = m\Y;
    end
    
    a = tuv(:,1)< 1 & tuv(:,1)>0;
    b = tuv(:,2)< 1 & tuv(:,2)>0;
    c = tuv(:,3)< 1 & tuv(:,3)>0;
    d = (tuv(:,2)+tuv(:,3))<=1;
    e = a&b&c&d;
    mul =min(tuv(e,1));
    
    if(~isempty(mul))
        outs(j,:) = Ia+(Ib-Ia)*mul;
    end
    
end
ind = sum(abs(outs),2)>0;
outs = outs(ind,:);


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
offset = S.offset;       
space = S.space;
wholehead = S.wholehead;
% Meshes
%--------------------------------------------------------------------------
scalp = gifti(D.inv{1}.mesh.tess_scalp);
cortex = gifti(D.inv{1}.mesh.tess_ctx);
Nf = size(scalp.faces,1);
lp = min(cortex.vertices(:,3));

% Face Positions(Downsampled)
%--------------------------------------------------------------------------

if(~wholehead)
    C= scalp.vertices(:,3)>(lp);
dscalp = spm_mesh_split(scalp,C);
else 
    dscalp =scalp;
end
p = [];
p.faces =dscalp.faces;
p.vertices = double(dscalp.vertices);
dscalp = reducepatch(p,.05);

[Nv,Nf]= spm_mesh_normals(dscalp,true);
posi = zeros(size(dscalp.faces));
cog = mean(scalp.vertices);

for i = 1:length(posi)
    whichVerts = dscalp.faces(i,:);
    cos = dscalp.vertices(whichVerts,:);
    posi(i,:) = mean(cos);   
end

% Face Orientations (Downsampled)
%--------------------------------------------------------------------------
ori = zeros(size(Nf));

for i = 1:length(Nf)
ad = posi(i,:) + 5*Nf(i,:);
subtrac = posi(i,:) - 5*Nf(i,:);

d1 = sum((cog - ad).^2);
d2 = sum((cog - subtrac).^2);

if(d2>d1)
    ori(i,:) = -Nf(i,:);
else
    ori(i,:) = Nf(i,:);
end
end

% Evenly sampled Positions on Downsampled scalp
%--------------------------------------------------------------------------
outs = poiGridSurface(dscalp,space,space*.8,space*1.25);

% Find what face we're on
%--------------------------------------------------------------------------
faceInd = zeros(size(outs,1),1);
for i = 1:length(outs)
    [mi,ind] = min(sqrt(sum(bsxfun(@minus,posi,outs(i,:)).^2,2)));
    faceInd(i)=ind;

end

% Project to used scalp
%--------------------------------------------------------------------------
pos = poiLineSurface(scalp,outs,ori(faceInd,:),space);


% Check if wholehead is requested
%--------------------------------------------------------------------------
if(~wholehead)
    C= scalp.vertices(:,3)>(lp);
    scalp = spm_mesh_split(scalp,C);
end
% Face Positions
%--------------------------------------------------------------------------
[Nv,Nf]= spm_mesh_normals(scalp,true);
posi = zeros(size(scalp.faces));
cog = mean(scalp.vertices);

for i = 1:length(posi)
    whichVerts = scalp.faces(i,:);
    cos = scalp.vertices(whichVerts,:);
    posi(i,:) = mean(cos);   
end

% Face Orientations 
%--------------------------------------------------------------------------
ori = zeros(size(Nf));

for i = 1:length(Nf)
ad = posi(i,:) + 5*Nf(i,:);
subtrac = posi(i,:) - 5*Nf(i,:);

d1 = sum((cog - ad).^2);
d2 = sum((cog - subtrac).^2);

if(d2>d1)
    ori(i,:) = -Nf(i,:);
else
    ori(i,:) = Nf(i,:);
end
end


% Find what face we're on
%--------------------------------------------------------------------------
faceInd = zeros(size(pos,1),1);
for i = 1:length(pos)
    [mi,ind] = min(sqrt(sum(bsxfun(@minus,posi,pos(i,:)).^2,2)));
    faceInd(i)=ind;
end


% add points if faces are empty and not near a sensor
%--------------------------------------------------------------------------

for i = 1:length(posi)
p=posi(i,:);
di= sort(sqrt(sum(bsxfun(@minus,pos,p).^2,2)));
if(di(1)>(space))
    pos = [pos;p];
    faceInd= [faceInd;i];
end
end
% return
%--------------------------------------------------------------------------

ori = ori(faceInd,:);
pos = pos+ori*offset;


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
 

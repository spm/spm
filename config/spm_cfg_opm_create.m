function create = spm_cfg_opm_create
% configuration file for creating OPM objects
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_cfg_opm_create.m 7429 2018-09-28 09:29:20Z spm $

%--------------------------------------------------------------------------
% Output Directory
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.val{1}  = {''};
outdir.help    = {'Files produced by this function will be written into this output directory. If no directory is given, images will be written to current working directory.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];
%--------------------------------------------------------------------------
% Output filename
%--------------------------------------------------------------------------
fname  = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'File name (.dat)';
fname.help    = {'File name, e.g. OPM.dat. If left empty this will be the name by default'};
fname.strtype = 's';
fname.num     = [1,Inf];
fname.val = {'OPM.dat'};

fifo = cfg_branch;
fifo.tag = 'fifo';
fifo.name = 'File management';
fifo.help = {'The output files and directories'};
fifo.val  = {outdir,fname};
%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
data        = cfg_files;
data.tag    = 'data';
data.name   = 'Labview .mat file';
data.filter = '.*.mat$';
data.num    = [1 1];
data.help   = {'Select the labview.mat file. If left empty an empty dataset will be created according to the simulation parameters'};
data.val = {{''}};
%--------------------------------------------------------------------------
% Sampling Frequency
%--------------------------------------------------------------------------
fs         = cfg_entry;
fs.tag     = 'fs';
fs.name    = 'Sampling Frequency';
fs.help    = {'The sampling frequency of the data in Hz'};
fs.strtype = 'r';
fs.num     = [1,1];
fs.val     = {1200};
%--------------------------------------------------------------------------
% Scale Factor
%--------------------------------------------------------------------------
scale        = cfg_entry;
scale.tag     = 'scale';
scale.name    = 'Scale Factor';
scale.help    = {'Scale factor (multiplied by data)to convert to units of fT'};
scale.strtype = 'r';
scale.num     = [1,1];
scale.val     = {1e6/2.7};

dataset = cfg_branch;
dataset.tag = 'dataset';
dataset.name = 'Dataset';
dataset.help = {'Provide basic information for dataset'};
dataset.val  = {data,fs,scale};

%--------------------------------------------------------------------------
% Pinout
%--------------------------------------------------------------------------
pinout = cfg_files;
pinout.tag = 'pinout';
pinout.name = 'Pinout File';
pinout.filter = '.*.txt';
pinout.num = [1 Inf];
pinout.help = {'Select the tab separated pinout .txt file.'};
pinout.val = {{''}};
%--------------------------------------------------------------------------
% OPM 2 Cast File
%--------------------------------------------------------------------------
sensorsUsed = cfg_files;
sensorsUsed.tag = 'sensorsUsed';
sensorsUsed.name = 'OPM to Cast File';
sensorsUsed.filter = '.*.txt';
sensorsUsed.num = [1 Inf];
sensorsUsed.help = {'Select the tab separated OPM to cast .txt file.(This is Scannercast and experiment specific)'};
sensorsUsed.val = {{''}};
%--------------------------------------------------------------------------
% Pos file
%--------------------------------------------------------------------------
pos = cfg_files;
pos.tag = 'pos';
pos.name = 'Scannercast specific position file ';
pos.filter = '.*.txt';
pos.num = [1 Inf];
pos.help = {'Select the tab separated position .txt file .'};
pos.val={{''}};

sens = cfg_branch;
sens.tag = 'sens';
sens.name = 'Sensor Information';
sens.help = {'Provide custom .txt files to describe sensor types and positions'};
sens.val  = {pinout,sensorsUsed,pos};

%--------------------------------------------------------------------------
% sMRI
%--------------------------------------------------------------------------
sMRI = cfg_files;
sMRI.tag = 'sMRI';
sMRI.name = 'Individual structural image';
sMRI.filter = 'image';
sMRI.ufilter = '.*';
sMRI.num     = [1 1];
sMRI.help = {'Select the subject''s structural image. Leave emptyp to create an OPM object without a forward model'};
sMRI.val = {{''}};
%--------------------------------------------------------------------------
% Meshes
%--------------------------------------------------------------------------
cortex = cfg_files;
cortex.tag = 'cortex';
cortex.name = 'Custom cortical mesh';
cortex.filter = 'mesh';
cortex.ufilter = '.*';
cortex.num     = [0 1];
cortex.help = {'Select the subject''s cortical mesh. Leave empty for default'};
cortex.val = {{''}};

iskull = cfg_files;
iskull.tag = 'iskull';
iskull.name = 'Custom inner skull mesh';
iskull.filter = 'mesh';
iskull.ufilter = '.*';
iskull.num     = [0 1];
iskull.help = {'Select the subject''s inner skull mesh. Leave empty for default'};
iskull.val = {{''}};

oskull = cfg_files;
oskull.tag = 'oskull';
oskull.name = 'Custom outer skull mesh';
oskull.filter = 'mesh';
oskull.ufilter = '.*';
oskull.num     = [0 1];
oskull.help = {'Select the subject''s outer skull mesh. Leave empty for default'};
oskull.val = {{''}};

scalp = cfg_files;
scalp.tag = 'scalp';
scalp.name = 'Custom scalp mesh';
scalp.filter = 'mesh';
scalp.ufilter = '.*';
scalp.num     = [0 1];
scalp.help = {'Select the subject''s scalp mesh. Leave empty for default'};
scalp.val = {{''}};

custom = cfg_branch;
custom.tag = 'custom';
custom.name = 'Custom meshes';
custom.help = {'Provide custom individual meshes as GIfTI files'};
custom.val  = {cortex, iskull, oskull, scalp};

meshres = cfg_menu;
meshres.tag = 'meshres';
meshres.name = 'Mesh resolution';
meshres.help = {'Specify the resolution of the cortical mesh'};
meshres.labels = {'coarse', 'normal', 'fine'};
meshres.values = {1, 2, 3};
meshres.val = {2};

meshing = cfg_branch;
meshing.tag = 'meshing';
meshing.name = 'Meshes';
meshing.help = {'Create head meshes for building the head model'};
meshing.val  = {sMRI,custom, meshres};


%--------------------------------------------------------------------------
% Volume Conducter
%--------------------------------------------------------------------------
voltype = cfg_menu;
voltype.tag = 'voltype';
voltype.name = 'MEG head model';
voltype.help = {'Select the head model type to use for MEG (if present)'};
voltype.labels = {'Single Sphere', 'MEG Local Spheres', 'Single Shell'};
voltype.values = {'Single Sphere', 'MEG Local Spheres', 'Single Shell'};
voltype.val = {'Single Shell'};


%--------------------------------------------------------------------------
% simulation parameters
%--------------------------------------------------------------------------

wholehead  = cfg_entry;
wholehead.tag     = 'wholehead';
wholehead.name    = 'Sensor Coverage';
wholehead.help    = {'If value is set to 1 then sensors will be placed all over scalp surface. '};
wholehead.strtype = 'r';
wholehead.num     = [1,1];
wholehead.val     = {1};

space  = cfg_entry;
space.tag     = 'space';
space.name    = 'Desired Space between Sensors';
space.help    = {'This is only an approximation, units are mm'};
space.strtype = 'r';
space.num     = [1,1];
space.val     = {25};

offset  = cfg_entry;
offset.tag     = 'offset';
offset.name    = 'Sensor Offset';
offset.help    = {'Scalp to sensor Distance'};
offset.strtype = 'r';
offset.num     = [1,1];
offset.val     = {6.5};

nSamples  = cfg_entry;
nSamples.tag     = 'nSamples';
nSamples.name    = 'Number of samples';
nSamples.help    = {''};
nSamples.strtype = 'r';
nSamples.num     = [1,1];
nSamples.val     = {1000};

lead  = cfg_entry;
lead.tag     = 'lead';
lead.name    = 'Compute Lead Field';
lead.help    = {'If value is set to 1 then a lead field will be compted and saved'};
lead.strtype = 'r';
lead.num     = [1,1];
lead.val     = {0};

simulation = cfg_branch;
simulation.tag = 'simulation';
simulation.name = 'Simulation Parameters';
simulation.help = {'Parameters for simulating OPM data. will be ignored if a dataset has alread been supplied'};
simulation.val  = {wholehead,space,offset,nSamples,lead};

%--------------------------------------------------------------------------
% simulation parameters
%--------------------------------------------------------------------------
create          = cfg_exbranch;
create.tag      = 'create';
create.name     = 'Create OPM object';
create.val      = {fifo,dataset,sens,meshing,voltype,simulation};
create.help     = {'Create/simulate OPM data. All arguments for this function are optional. It allows for either a conversion of raw data to valid MEG object or for the simulation of OPM data on a template brain or an individual brain with customisable sensor configurations.'}';
create.prog     = @opm_create;
create.vout     = @vout_opm_create;
create.modality = {'EEG'};


%==========================================================================
function out = opm_create(job)
% construct the S struct

% datset parameters
S=[];
S.data=job.dataset.data{1};
S.fs=job.dataset.fs;
S.scale= job.dataset.scale;

% sensor information
S.pinout =job.sens.pinout{1};
S.pos=job.sens.pos{1};
S.sensorsUsed= job.sens.sensorsUsed{1};

%  meshses
S.cortex = job.meshing.custom.cortex;
S.iskull = job.meshing.custom.iskull;
S.oskull = job.meshing.custom.oskull;
S.scalp  = job.meshing.custom.scalp;
S.meshres = job.meshing.meshres;
S.sMRI = job.meshing.sMRI{1};
S.voltype = job.voltype;

%  simulation
S.lead= job.simulation.lead;
S.nSamples = job.simulation.nSamples;
S.offset = job.simulation.offset;
S.space = job.simulation.space;
S.wholehead = job.simulation.wholehead;

% set filename(will be udpate later if necessary)
S.fname = fullfile(job.fifo.outdir{1},job.fifo.fname);

% check if field should be removed
argumentFields = fields(S);
numFields = length(argumentFields);
notUsed = zeros(numFields,1);

for i = 1:numFields
    notUsed(i) = iscell(S.(argumentFields{i}));
end

% remove the unused fields
fieldsToRemove = {argumentFields{boolean(notUsed),1}};
S = rmfield(S,fieldsToRemove);


% check if dataset has been provided, load and rename. 
if(isfield(S, 'data'))
    load(job.dataset.data{1});
    S.data= data.B';
    trigs = (size(data.decimalTrigs,2)+size(data.binaryTrigs,2))>0;
    if(trigs)
        S.trig= [data.decimalTrigs,data.binaryTrigs]';
    end
    [a,b,c]=fileparts(job.dataset.data{1});
    outfile= fullfile(a,['SPM_',b,'.dat']);
    S.fname = outfile;
end

% run the main function 
out.D= spm_opm_create(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};


%==========================================================================
function dep = vout_opm_create(job)
% return dependencies
dep = cfg_dep;
dep.sname = 'OPM Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'OPM Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});

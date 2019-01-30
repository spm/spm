function create = spm_cfg_opm_create
% configuration file for creating OPM objects
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_cfg_opm_create.m 7521 2019-01-30 18:16:03Z tim $

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
data        = cfg_files;
data.tag    = 'data';
data.name   = 'OPM dataset';
data.filter = '.*.bin$';
data.num    = [1 1];
data.help   = {'Select the meg.bin file. If left empty an empty dataset will be created according to the simulation parameters'};
data.val = {{''}};

chans = cfg_files;
chans.tag = 'chan';
chans.name = 'channels.tsv file';
chans.filter = '.*.tsv';
chans.num = [1 Inf];
chans.help = {'Select the channels .tsv file. The format of this file should conform to the BIDS standard for channels.tsv files'};
chans.val = {{''}};

meg = cfg_files;
meg.tag = 'meg';
meg.name = 'meg.json file';
meg.filter = '.*.json';
meg.num = [1 Inf];
meg.help = {'Select the meg  json file.  The format of this file should conform to the BIDS standard for meg.json files'};
meg.val = {{''}};
%--------------------------------------------------------------------------
% Sensor Branch
%--------------------------------------------------------------------------
sens = cfg_branch;
sens.tag = 'sens';
sens.name = 'Sensor Info';
sens.help = {'Input arguements required for sensor level analysis'};
sens.val  = {data,chans,meg};

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
offset.help    = {'Scalp to sensor Distance, units are mm'};
offset.strtype = 'r';
offset.num     = [1,1];
offset.val     = {6.5};

nSamples  = cfg_entry;
nSamples.tag     = 'nSamples';
nSamples.name    = 'Number of samples';
nSamples.help    = {''};
nSamples.strtype = 'r';
nSamples.num     = [1,1];
nSamples.val     = {1};

%--------------------------------------------------------------------------
% simulation branch
%--------------------------------------------------------------------------
simulation = cfg_branch;
simulation.tag = 'simulation';
simulation.name = 'Simulation Parameters';
simulation.help = {'Parameters for simulating OPM data. will be ignored if a dataset has alread been supplied'};
simulation.val  = {wholehead,space,offset,nSamples};


%--------------------------------------------------------------------------
% Source paramters
%--------------------------------------------------------------------------
positions = cfg_files;
positions.tag = 'pos';
positions.name = 'Sensor positins';
positions.filter = '.*.tsv';
positions.num = [1 Inf];
positions.help = {'tsv file giving coordinates for labelled sensors'};
positions.val={{''}};

coordsystem = cfg_files;
coordsystem.tag = 'coord';
coordsystem.name = 'Coordinate systems';
coordsystem.filter = '.*.json';
coordsystem.num = [1 Inf];
coordsystem.help = {'json file describing hte coordinate systesm of fiducials and sensors.The format of this file should conform to the BIDS standard for coordsystem.json files'};
coordsystem.val={{''}};

sMRI = cfg_files;
sMRI.tag = 'sMRI';
sMRI.name = 'Individual structural image';
sMRI.filter = 'image';
sMRI.ufilter = '.*';
sMRI.num     = [1 1];
sMRI.help = {'Select the subject''s structural image.'};
sMRI.val = {{''}};

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

voltype = cfg_menu;
voltype.tag = 'voltype';
voltype.name = 'MEG head model';
voltype.help = {'Select the head model type to use for MEG (if present)'};
voltype.labels = {'Single Sphere', 'MEG Local Spheres', 'Single Shell'};
voltype.values = {'Single Sphere', 'MEG Local Spheres', 'Single Shell'};
voltype.val = {'Single Shell'};

lead  = cfg_entry;
lead.tag     = 'lead';
lead.name    = 'Compute Lead Field';
lead.help    = {'If value is set to 1 then a lead field will be compted and saved'};
lead.strtype = 'r';
lead.num     = [1,1];
lead.val     = {0};

%--------------------------------------------------------------------------
% source branch
%--------------------------------------------------------------------------
source = cfg_branch;
source.tag = 'source';
source.name = 'Source level info';
source.help = {'Input arguments necessary for performing source space analysis'};
source.val  = {coordsystem,positions,sMRI,meshres,custom,voltype,lead};

%--------------------------------------------------------------------------
% Create branch
%--------------------------------------------------------------------------
create          = cfg_exbranch;
create.tag      = 'create';
create.name     = 'Create OPM object';
create.val      = {sens,simulation,source};
create.help     = {'Create/simulate OPM data. All arguments for this function are optional. It allows for either a conversion of raw data to valid MEG object or for the simulation of OPM data on a template brain or an individual brain with customisable sensor configurations.'}';
create.prog     = @opm_create;
create.vout     = @vout_opm_create;
create.modality = {'EEG'};


%==========================================================================
function out = opm_create(job)
% construct the S struct

% datset parameters
S=[];
S.data=job.sens.data{1};
S.meg=job.sens.meg;
S.channels =job.sens.chan{1};

%  meshses
S.cortex = job.source.custom.cortex{1};
S.iskull = job.source.custom.iskull{1};
S.oskull = job.source.custom.oskull{1};
S.scalp  = job.source.custom.scalp{1};
S.meshres = job.source.meshres;
S.sMRI = job.source.sMRI{1};
S.voltype = job.source.voltype;
S.lead= job.source.lead;
S.coordsystem= job.source.coord{1};
S.positions=job.source.pos{1};

%  simulation
S.nSamples = job.simulation.nSamples;
S.offset = job.simulation.offset;
S.space = job.simulation.space;
S.wholehead = job.simulation.wholehead;

% check if field should be removed
argumentFields = fields(S);
numFields = length(argumentFields);
notUsed = zeros(numFields,1);

for i = 1:numFields
    notUsed(i) = strcmp(S.(argumentFields{i}),'');
end

% remove the unused fields
fieldsToRemove = {argumentFields{boolean(notUsed),1}};
S = rmfield(S,fieldsToRemove);

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

function extract = spm_cfg_eeg_inv_extract
% configuration file for extracting source data from imaging source
% reconstruction
%_______________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_inv_extract.m 5377 2013-04-02 17:07:57Z vladimir $

D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG datasets';
D.filter = 'mat';
D.num = [1 Inf];
D.help = {'Select the M/EEG mat files.'};

val = cfg_entry;
val.tag = 'val';
val.name = 'Inversion index';
val.strtype = 'n';
val.help = {'Index of the cell in D.inv where the inversion results are be stored.'};
val.val = {1};

xyz = cfg_entry;
xyz.tag = 'xyz';
xyz.name = 'Source location';
xyz.strtype = 'r';
xyz.num = [1 3];
xyz.help = {'Source location (in MNI coordinates)'};

label = cfg_entry;
label.tag = 'label';
label.name = 'Source label';
label.strtype = 's';
label.num = [1 Inf];
label.help = {'Label for the source channel in the output file'};

source = cfg_branch;
source.tag = 'source';
source.name = 'Source';
source.val = {label, xyz};

sources = cfg_repeat;
sources.tag = 'sources';
sources.name = 'Sources to extract';
sources.values  = {source};
sources.num     = [1 Inf];
sources.help = {'Specify sources to extract data from'};

rad = cfg_entry;
rad.tag = 'rad';
rad.name = 'VOI radius';
rad.strtype = 'r';
rad.num = [1 1];
rad.val = {5};
woi.help = {'Radius around each location to extract an eigenvariate from (mm).'};

type = cfg_menu;
type.tag = 'type';
type.name = 'What to extract';
type.help = {'What to extract: evoked activity or single trials.'};
type.labels = {'Evoked', 'Single trials'};
type.values = {'evoked', 'trials'};
type.val = {'trials'};

fname = cfg_entry;
fname.tag = 'fname';
fname.name = 'Output dataset name';
fname.strtype = 's';
fname.num = [0 Inf];
fname.val = {''};
fname.help = {'Output file name (empty for default)'};

extract = cfg_exbranch;
extract.tag = 'extract';
extract.name = 'Source extraction';
extract.val = {D, val, sources, rad, type, fname};
extract.help = {'Extract source data from the results of inverse source reconstruction'};
extract.prog = @run_extract;
extract.vout = @vout_extract;
extract.modality = {'EEG'};

function  out = run_extract(job)

source       = [];
source.XYZ   = cat(1, job.source.xyz);
source.label = {job.source.label};
source.rad   = job.rad;
source.type  = job.type;

if ~isempty(job.fname)
    source.fname = job.fname;
end

out.D = {};

for i = 1:numel(job.D)
    D = spm_eeg_load(job.D{i});
    
    D.val = job.val;      
            
    if ~isfield(D.inv{D.val}, 'inverse') || ~isfield(D.inv{D.val}.inverse, 'J')
        error('Imaging source reconstruction is missing for subject %d.', i);
    end       
    
    D.inv{D.val}.source = source;
    
    Ds = spm_eeg_inv_extract(D);
    
    out.D{i, 1} = fullfile(Ds.path, Ds.fname);      
end

function dep = vout_extract(job)
% Output is always in field "D", no matter how job is structured
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) extracted source data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

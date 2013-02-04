function simulate = spm_cfg_eeg_inv_simulate
% configuration file for configuring imaging source inversion
% reconstruction
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_inv_simulate.m 5230 2013-02-04 14:03:21Z gareth $

D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG datasets';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the template M/EEG mat file.'};

prefix = cfg_entry;
prefix.tag = 'prefix';
prefix.name = 'Output file prefix';
prefix.strtype = 's';
prefix.val = {'sim_'};

val = cfg_entry;
val.tag = 'val';
val.name = 'Inversion index';
val.strtype = 'n';
val.help = {'Index of the cell in D.inv where the forward model can be found.'};
val.val = {1};

all = cfg_const;
all.tag = 'all';
all.name = 'All';
all.val  = {1};

condlabel = cfg_entry;
condlabel.tag = 'condlabel';
condlabel.name = 'Condition label';
condlabel.strtype = 's';
condlabel.val = {'faces'};

conditions = cfg_repeat;
conditions.tag = 'conditions';
conditions.name = 'Conditions';
conditions.help = {'Specify the labels of the conditions to be replaced by simulated data'};
conditions.num  = [1 Inf];
conditions.values  = {condlabel};
conditions.val = {condlabel};

whatconditions = cfg_choice;
whatconditions.tag = 'whatconditions';
whatconditions.name = 'What conditions to include?';
whatconditions.values = {all, conditions};
whatconditions.val = {all};

% standard = cfg_const;
% standard.tag = 'standard';
% standard.name = 'Standard';
% standard.help = {'Use default settings for the simulation'};
% standard.val  = {1};
% %
% invtype = cfg_menu;
% invtype.tag = 'invtype';
% invtype.name = 'Inversion type';
% invtype.help = {'Select the desired inversion type'};
% invtype.labels = {'GS', 'ARD', 'MSP (GS+ARD)' 'COH', 'IID'};
% invtype.values = {'GS', 'ARD', 'MSP', 'LOR', 'IID'};
% invtype.val = {'GS'};

woi = cfg_entry;
woi.tag = 'woi';
woi.name = 'Time window';
woi.strtype = 'r';
woi.num = [1 2];
woi.val = {[0.1 0.4]};
woi.help = {'Time window in which to simulate data (s)'};


foi = cfg_entry;
foi.tag = 'foi';
foi.name = 'Frequencies (Hz) ';
foi.strtype = 'r';
foi.num = [Inf 1];
foi.val = {[10 ; 20]};
foi.help = {'Enter frequencies of sources in Hz'};


dipmom = cfg_entry;
dipmom.tag = 'dipmom';
dipmom.name = 'Dipole moments (nAm) ';
dipmom.strtype = 'r';
dipmom.num = [Inf 1];
dipmom.val = {[20;20]};
dipmom.help = {'Enter dipole moments for sources (nAm)'};

%
% nsources = cfg_entry;
% nsources.tag = 'nsources';
% nsources.name = 'Number of sources';
% nsources.strtype = 'i';
% nsources.num = [1 1];
% nsources.val = {[2]};
% nsources.help = {'Number of simulated patches'};
%
% patchind = cfg_entry;
% patchind.tag = 'patchind';
% patchind.name = 'Number of iterations';
% patchind.strtype = 'i';
% patchind.num = [1 2];
% patchind.val = {[100 150]};
% patchind.help = {'Number of times the inversion will be run using a random set of patches each time'};
%
% nsmodes = cfg_entry;
% nsmodes.tag = 'nsmodes';
% nsmodes.name = 'Number of spatial modes';
% nsmodes.strtype = 'i';
% nsmodes.num = [1 1];
% nsmodes.val = {[100]};
% nsmodes.help = {'Number of spatial modes'};
%
% ntmodes = cfg_entry;
% ntmodes.tag = 'ntmodes';
% ntmodes.name = 'Number of temporal modes';
% ntmodes.strtype = 'i';
% ntmodes.num = [1 1];
% ntmodes.val = {[4]};
% ntmodes.help = {'Number of temporal modes'};

%
% priorsmask  = cfg_files;
% priorsmask.tag = 'priorsmask';
% priorsmask.name = 'Priors file';
% priorsmask.filter = '(.*\.gii$)|(.*\.mat$)|(.*\.nii(,\d+)?$)|(.*\.img(,\d+)?$)';
% priorsmask.num = [0 1];
% priorsmask.help = {'Select a mask or a mat file with priors.'};
% priorsmask.val = {{''}};
%
% space = cfg_menu;
% space.tag = 'space';
% space.name = 'Prior image space';
% space.help = {'Space of the mask image.'};
% space.labels = {'MNI', 'Native'};
% space.values = {1, 0};
% space.val = {1};
%
% priors = cfg_branch;
% priors.tag = 'priors';
% priors.name = 'Source priors';
% priors.help = {'Restrict solutions to pre-specified VOIs'};
% priors.val  = {priorsmask, space};

locs  = cfg_entry;
locs.tag = 'locs';
locs.name = 'Source locations (or patch indices)';
locs.strtype = 'r';
locs.num = [Inf 3];
locs.help = {'Input mni source locations (mm) as n x 3 matrix, or patch indices only in 1st column'};
locs.val = {[ 53.7203  -25.7363    9.3949;  -52.8484  -25.7363    9.3949]};
%
% radius = cfg_entry;
% radius.tag = 'radius';
% radius.name = 'Radius of VOI (mm)';
% radius.strtype = 'r';
% radius.num = [1 1];
% radius.val = {32};
%
% restrict = cfg_branch;
% restrict.tag = 'restrict';
% restrict.name = 'Restrict solutions';
% restrict.help = {'Restrict solutions to pre-specified VOIs'};
% restrict.val  = {locs, radius};
%
% custom = cfg_branch;
% custom.tag = 'custom';
% custom.name = 'Custom';
% custom.help = {'Define custom settings for the inversion'};
% custom.val  = {invtype, woi, foi, hanning,nsources,patchind,nsmodes,ntmodes, priors, restrict};

% isstandard = cfg_choice;
% isstandard.tag = 'isstandard';
% isstandard.name = 'Inversion parameters';
% isstandard.help = {'Choose whether to use standard or custom inversion parameters.'};
% isstandard.values = {standard, custom};
% isstandard.val = {standard};
%
% modality = cfg_menu;
% modality.tag = 'modality';
% modality.name = 'Select modalities';
% modality.help = {'Select modalities for the inversion (only relevant for multimodal datasets).'};
% modality.labels = {'All', 'EEG', 'MEG', 'MEGPLANAR', 'EEG+MEG', 'MEG+MEGPLANAR', 'EEG+MEGPLANAR'};
% modality.values = {
%     {'All'}
%     {'EEG'}
%     {'MEG'}
%     {'MEGPLANAR'}
%     {'EEG', 'MEG'}
%     {'MEG', 'MEGPLANAR'}
%     {'EEG', 'MEGPLANAR'}
%     }';
% modality.val = {{'All'}};

simulate = cfg_exbranch;
simulate.tag = 'simulate';
simulate.name = 'simulation of sources on cortex';
simulate.val = {D, val, prefix, condlabel, whatconditions,locs,foi,woi,dipmom};
simulate.help = {'Run simulation'};
simulate.prog = @run_simulation;
simulate.vout = @vout_simulation;
simulate.modality = {'EEG'};

function  out = run_simulation(job)


D = spm_eeg_load(job.D{1});

sim = [];
if isfield(job.whatconditions, 'condlabel')
    sim.trials = job.whatconditions.condlabel;
end
if numel(job.D)>1,
    error('Simulattion routine only meant to replace data for single subjects');
end;




sim.woi  = [max(job.woi(1), D.time(1)) min(job.woi(2), D.time(end))];

Nsources=size(job.foi,1);
if size(job.locs,1)~=Nsources,
    error('Number of locations must equal number of frequencies specified');
end;







[mod, list] = modality(D, 1, 1);
if ~strcmp(mod, 'MEG')
    error('only suitable for pure MEG data at the moment');
end;

D = {};



for i = 1:numel(job.D) %% only set up for one subject at the moment but leaving this for the future
    D{i} = spm_eeg_load(job.D{i});
    D{i}.val = job.val;
    
    D{i}.con = 1;
    if ~isfield(D{i}, 'inv')
        error(sprintf('Forward model is missing for subject %d', i));
    elseif  numel(D{i}.inv)<D{i}.val || ~isfield(D{i}.inv{D{i}.val}, 'forward')
        if D{i}.val>1 && isfield(D{i}.inv{D{i}.val-1}, 'forward')
            D{i}.inv{D{i}.val} = D{i}.inv{D{i}.val-1};
            warning(sprintf('Duplicating the last forward model for subject %d', i));
        else
            error(sprintf('Forward model is missing for subject %d', i));
        end
    end
    
end; % for i


mnimesh=[]; %% using mesh defined in forward model at the moment
SmthInit=[]; %% leave patch size as default for now
keep_leads=[]; %% recompute lead fields each time

SIdipmom=job.dipmom*1e-9; %% put dipole moment into Am rather than nAm
[D,meshsourceind,signal]=spm_eeg_simulate(D,job.prefix, job.locs,job.foi,job.woi,SIdipmom,mnimesh,keep_leads,SmthInit);

if ~iscell(D)
    D = {D};
end

for i = 1:numel(D)
    save(D{i});
end

out.D = job.D;

function dep = vout_simulation(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after simulation';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});


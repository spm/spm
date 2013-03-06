function simulate = spm_cfg_eeg_inv_simulate
% configuration file for configuring imaging source inversion
% reconstruction
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_inv_simulate.m 5307 2013-03-06 17:15:45Z gareth $

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


woi = cfg_entry;
woi.tag = 'woi';
woi.name = 'Time window';
woi.strtype = 'r';
woi.num = [1 2];
woi.val = {[100  400]};
woi.help = {'Time window in which to simulate data (ms)'};


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


setSNR = cfg_entry;
setSNR.tag = 'setSNR';
setSNR.name = 'Sensor level SNR (dBs)';
setSNR.strtype = 'r';
setSNR.num = [1 1];
setSNR.val = {[0]};
setSNR.help = {'Enter sensor level SNR=20*log10(rms source/ rms noise)'};


isSNR = cfg_choice;
isSNR.tag = 'isSNR';
isSNR.name = 'Set SNR or set moments';
isSNR.help = {'Choose whether to a fixed SNR or specify source moments in nAm'};
isSNR.values = {setSNR, dipmom};
isSNR.val = {setSNR};

locs  = cfg_entry;
locs.tag = 'locs';
locs.name = 'Source locations (or patch indices)';
locs.strtype = 'r';
locs.num = [Inf 3];
locs.help = {'Input mni source locations (mm) as n x 3 matrix, or patch indices only in 1st column'};
locs.val = {[ 53.7203  -25.7363    9.3949;  -52.8484  -25.7363    9.3949]};

simulate = cfg_exbranch;
simulate.tag = 'simulate';
simulate.name = 'simulation of sources on cortex';
simulate.val = {D, val, prefix, whatconditions,locs,foi,woi,isSNR};
simulate.help = {'Run simulation'};
simulate.prog = @run_simulation;
simulate.vout = @vout_simulation;
simulate.modality = {'EEG'};

function  out = run_simulation(job)


D = spm_eeg_load(job.D{1});


trialind=[];
if isfield(job.whatconditions, 'condlabel')
    trialind =D.indtrial( job.whatconditions.condlabel);
    if isempty(trialind),
        error('No such condition found');
    end; 
end
if numel(job.D)>1,
    error('Simulation routine only meant to replace data for single subjects');
end;



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


if isfield(job.isSNR,'dipmom'),
    SIdipmom=job.isSNR.dipmom*1e-9; %% put dipole moment into Am rather than nAm
    SNRdB=[];
    whitenoise=[]; %% need to add this later
else
    SNRdB=job.isSNR.setSNR;
    SIdipmom=[];
    whitenoise=[];
end;

[D,meshsourceind,signal]=spm_eeg_simulate(D,job.prefix, job.locs,job.foi,job.woi./1000,SIdipmom,whitenoise,SNRdB,trialind,mnimesh,SmthInit);

if ~iscell(D)
    D = {D};
end

for i = 1:numel(D)
    save(D{i});
end

out.D = {D{1}.fname};

function dep = vout_simulation(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after simulation';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});


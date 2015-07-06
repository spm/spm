function simulate = spm_cfg_eeg_inv_simulate
% configuration file for configuring imaging source inversion
% reconstruction
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_inv_simulate.m 6494 2015-07-06 10:23:04Z gareth $

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
foi.name = 'Frequencies of sinusoid (Hz) ';
foi.strtype = 'r';
foi.num = [Inf 1];
foi.val = {[10 ; 20]};
foi.help = {'Enter frequencies of sources in Hz'};

fband = cfg_entry;
fband.tag = 'fband';
fband.name = 'Bandwidth of Gaussian orthogonal white signals in (Hz) ';
fband.strtype = 'r';
fband.num = [ 1 2];
fband.val = {[10  40]};
fband.help = {'Enter frequencies over which random signal exists in Hz'};


fromfile = cfg_files;
fromfile.tag = 'fromfile';
fromfile.name = 'Variable s in a .mat file';
fromfile.filter = 'mat';
fromfile.num = [1 1];
fromfile.help = {'Select the .mat file containing the variable s'};


isSin = cfg_choice;
isSin.tag = 'isSin';
isSin.name = 'How to generate signals ?';
isSin.help = {'Choose whether to simulate a sinusoidal signal OR a random (gaussian white) filtered (and orthogonal) signals OR to load time series from a matfile (s(n,:) is nth dipole time series (all timeserie will be normalized to unit amplitude, and scaled only by dipole moment)'};
isSin.values = {foi, fband,fromfile};
isSin.val = {foi};


dipmom = cfg_entry;
dipmom.tag = 'dipmom';
dipmom.name = 'Dipole moment  ';
dipmom.strtype = 'r';
dipmom.num = [Inf Inf];
dipmom.val = {[10 10;10 5]};
dipmom.help = {' EITHER (for perfect dipoles) enter dipole moment in nAm in x,y,z as nAm OR (for surfaces) enter total dipole moment (again nAm), followed by spatial extent (FWHM of a Gaussian) in mm. Note only relative value will be used if SNR is specified later'};

whitenoise = cfg_entry;
whitenoise.tag = 'whitenoise';
whitenoise.name = 'White noise in rms femto Tesla ';
whitenoise.strtype = 'r';
whitenoise.num = [1 1];
whitenoise.val = {[100]};
whitenoise.help = {'White noise in the recording bandwidth (rms fT)'};


setSNR = cfg_entry;
setSNR.tag = 'setSNR';
setSNR.name = 'Sensor level SNR (dBs)';
setSNR.strtype = 'r';
setSNR.num = [1 1];
setSNR.val = {[0]};
setSNR.help = {'Enter sensor level SNR=20*log10(rms source/ rms noise)'};

locs  = cfg_entry;
locs.tag = 'locs';
locs.name = 'Source locations (mm)';
locs.strtype = 'r';
locs.num = [Inf 3];
locs.help = {'Input mni source locations (mm) as n x 3 matrix'};
locs.val = {[ 53.7203  -25.7363    9.3949;  -52.8484  -25.7363    9.3949]};


setsources = cfg_branch;
setsources.tag = 'setsources';
setsources.name = 'Set sources';
setsources.help = {'Define sources and locations (in MNI space)'};
setsources.val  = {woi, isSin, dipmom, locs};

setinversion = cfg_const;
setinversion.tag = 'isinversion';
setinversion.name = 'Use current density estimate';
setinversion.help = {'Use the current density estimate at the inversion index to produce MEG data'};
setinversion.val  = {1};

isinversion = cfg_choice;
isinversion.tag = 'isinversion';
isinversion.name = 'Use inversion or define sources';
isinversion.help = {'Use existing current density estimate to generate data or generate from defined sources'};
isinversion.values = {setinversion, setsources};
isinversion.val = {setsources};



isSNR = cfg_choice;
isSNR.tag = 'isSNR';
isSNR.name = 'Set SNR or set white noise level';
isSNR.help = {'Choose whether to a fixed SNR or specify system noise level'};
isSNR.values = {setSNR, whitenoise};
isSNR.val = {setSNR};


simulate = cfg_exbranch;
simulate.tag = 'simulate';
simulate.name = 'Simulation of sources';
simulate.val = {D, val, prefix, whatconditions,isinversion,isSNR};
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



if size(modality(D),1)>1,
    error('only suitable for single modality data at the moment');
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


if isfield(job.isSNR,'whitenoise'),
    whitenoisefT=job.isSNR.whitenoise; %% internal units as femto Tesla
    SNRdB=[];
else
    SNRdB=job.isSNR.setSNR;
    whitenoisefT=[];
end;

if isfield(job.isinversion,'setsources'), %% defining individual sources
    
    %%%
    Nsources=size(job.isinversion.setsources.locs,1)
    
    if (size(job.isinversion.setsources.dipmom,1)~=Nsources),
        error('Number of locations must equal number of moments specified');
    end;
    
    mnimesh=[]; %% using mesh defined in forward model at the moment
    
    
    ormni=[]; %% dipoles will get orientations from cortical mesh
    dipmom=job.isinversion.setsources.dipmom;
    if size(job.isinversion.setsources.dipmom,2)==3, %% dipole orientation is specified
        disp('Simulating dipoles without reference to cortical surface');
        for j=1:size(dipmom,1),
            ormni(j,:)=dipmom(j,:)./sqrt(dot(dipmom(j,:),dipmom(j,:))); %% unit orientation in MNI space
            nAmdipmom(j)=sqrt(dot(dipmom(j,:),dipmom(j,:))); % magnitude of dipole
        end;
    else %% only one moment parameter given
        nAmdipmom=dipmom(:,1); %% total momnent in nAm
        dipfwhm=[];
        if size(dipmom,1)==2,
            dipfwhm=dipmom(:,2); %% fhwm in mm
        end;
    end;
    
    woi=job.isinversion.setsources.woi./1000;
    timeind = intersect(find(D{1}.time>woi(1)),find(D{1}.time<=woi(2)));
    simsignal = zeros(Nsources,length(timeind));
    
    if isfield(job.isinversion.setsources.isSin,'fromfile'),
        %% Simulate orthogonal Gaussian signals
        filename=cell2mat(job.isinversion.setsources.isSin.fromfile);
        a=load(filename,'s');
        
        if size(a.s,1)~=Nsources,
            error('size of simulated data in file does not match number of sources');
        end;
        
        if size(a.s,2)~=length(timeind),
            warning(sprintf('sim signal from file is not same length as time window (%d vs %d samples) truncating or padding with zeros',size(a.s,2),length(timeind)));
        end;
        usesamples=1:min(length(timeind),size(a.s,2));
        simsignal(:,usesamples)=a.s(:,usesamples);
        
    end; % if isfield fband
    
    if isfield(job.isinversion.setsources.isSin,'fband'),
        %% Simulate orthogonal Gaussian signals
        
        simsignal=randn(Nsources,length(timeind));
        %% filter to bandwidth
        simsignal=ft_preproc_lowpassfilter(simsignal,D{1}.fsample,job.isinversion.setsources.isSin.fband(2),2);
        simsignal=ft_preproc_highpassfilter(simsignal,D{1}.fsample,job.isinversion.setsources.isSin.fband(1),2);
        [u,s,v]=svd(simsignal*simsignal');
        simsignal=u*simsignal; %% orthogonalise all signals
    end; % if isfield fband
    
    %  simsignal=simsignal.*repmat(nAmdipmom,1,size(simsignal,2)); %% now scale by moment
    
    if isfield(job.isinversion.setsources.isSin,'foi'),
        %% simulate sinusoids
        sinfreq=job.isinversion.setsources.isSin.foi;
        % Create the waveform for each source
        
        for j=1:Nsources                % For each source
            simsignal(j,:)=sin((D{1}.time(timeind)- D{1}.time(min(timeind)))*sinfreq(j)*2*pi);
        end; % for j
        
        
    end; %% if isfield foi
    
    simsignal=simsignal./repmat(std(simsignal'),size(simsignal,2),1)'; %% Set sim signal to have unit variance
    %[D,meshsourceind,signal]=spm_eeg_simulate(D,job.prefix, job.isinversion.setsources.locs,simsignal,woi,whitenoisefT,SNRdB,trialind,mnimesh,SmthInit);
    figure;
    plot(D{1}.time(timeind),simsignal);
    xlabel('time');
    ylabel('Normalized amplitude (to be scaled by dip moment later)');
    legend(num2str([1:Nsources]'));
    [D,meshsourceind]=spm_eeg_simulate(D,job.prefix,job.isinversion.setsources.locs,simsignal,ormni,woi,whitenoisefT,SNRdB,trialind,mnimesh,dipfwhm,nAmdipmom);
    
else %% simulate sources based on inversion
    if ~isfield(D{i}.inv{job.val},'inverse'),
        error('There is no solution defined for these data at that index');
    end;
    
    [D]=spm_eeg_simulate_frominv(D,job.prefix,job.val,whitenoisefT,SNRdB,trialind);
end;



if ~iscell(D)
    D = {D};
end

for i = 1:numel(D)
    save(D{i});
end


fullname=[D{1}.path filesep D{1}.fname];
out.D = {fullname};

function dep = vout_simulation(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after simulation';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});


function S = spm_config_meeg
% configuration file for MEEG preprocessing
%_______________________________________________________________________
% Copyright (C) 2006 Wellcome Department of Neuroimaging

%% ------------ Files ------------
dataset = struct('type','files','name','Select the dataset file','tag','dataset',...
    'filter','*','num',1,...
    'help',{{'Select the data file'}});

% ------------Trials---------------

eventtype = struct('type','entry','name','Event type','tag','eventtype',...
    'strtype', 's', 'num', [1 1], 'val', {{'gui'}}, 'help', {{'Type of the event for trial selection'}});

eventvalue = struct('type','entry','name','Event value','tag','eventvalue',...
    'strtype', 's', 'num', [1 1], 'val',  {{'gui'}},'help', {{'Value of the event for trial selection'}});

prestim = struct('type','entry','name','Pre-trigger time','tag','prestim',...
    'strtype', 'r', 'num', [1 1], 'help', {{'Pre-trigger latency in seconds'}});

poststim = struct('type','entry','name','Post-trigger time','tag','poststim',...
    'strtype', 'r', 'num', [1 1], 'help', {{'Post-trigger latency in seconds'}});

ignorelim=struct('type','menu','name','Ignore trial limits in the data','tag','ignorelim',...
    'labels', {{'Yes', 'No'}}, 'values', {{'yes', 'no'}}, 'val', {{'no'}});

ignorelim.help= [...
    {'In some file formats there can be events that already '},...
    {'contain information about the trial limits. If such event'},...
    {'is chosen SPM will disregard the user defined pre- and post-'},...
    {'stimulus times. Choose ''yes'' in this menu to override this'},...
    {'default behavior and ignore the information from the file.'},...
    {'Leave it as ''no'' unless you know well what you are doing.'}];


savetrl = struct('type','menu','name','Save TRL to file','tag','savetrl',...
    'labels', {{'Yes', 'No'}}, 'values', {{'yes', 'no'}}, 'val', {{'no'}},...
    'help', {{'Choose whether to save the TRL matrix to a mat or ASCII file'}});


sub_trialdef = struct('type','branch','name','Define based on events','tag','sub_trialdef',...
    'val',{{eventtype, eventvalue, prestim, poststim, ignorelim, savetrl}});


loadcont=struct('type','const','name','Load continuous data','tag','loadcont', 'val', {{' '}},...
    'help', {{'Load continuous data without epoching.'}});


trlfile = struct('type','files','name','Define by TRL file','tag','trlfile',...
    'filter','*', 'num',1, 'help',{{'Select an ASCII or matlab file with TRL structure (see definitrial() help)'}});


trials = struct('type','choice','name','Trial definitions','tag','trials',...
    'values',{{sub_trialdef, trlfile, loadcont}}, 'val', {{sub_trialdef}});

trials.help=[...
    {'Choose whether to define trials semi-automatically based on events in the file'}, ...
    {'or using iser defined trl matrix (see definitrial() help)'}];

forcecont = struct('type','menu','name','Force treating data as continuous','tag','forcecont',...
    'labels', {{'Yes', 'No'}}, 'values', {{'yes', 'no'}}, 'val', {{'no'}},...
    'help', {{'Force preprocessing to read across trial borders.'}});

%% ------------ Channels ------------

channel= struct('type','entry','name','Select channels','tag','channel',...
    'strtype', 's', 'val', {{'all'}}, 'num', [1 1]);

channel.help = [...
    {'''gui''     a graphical user interface will pop up to select the channels'},...
    {'''all''     is replaced by all channels in the datafile'},...
    {'''MEG''     is replaced by all channels in the CTF datafile starting with ''M'''},...
    {'''EEG''     is replaced by all channels in the CTF datafile starting with  ''EEG'''},...
    {'''EEG1020'' is replaced by ''Fp1'', ''Fpz'', ''Fp2'', ''F7'', ''F3'','} ,...
    {'''EOG''     is replaced by all recognized EOG channels'},...
    {'''EMG''     is replaced by all channels in the datafile starting with ''EMG'''},...
    {'''lfp''     is replaced by all channels in the datafile starting with ''lfp'''},...
    {'''mua''     is replaced by all channels in the datafile starting with ''mua'''},...
    {'''spike''   is replaced by all channels in the datafile starting with ''spike'''},...
    {' 10      is replaced by the 10th channel in the datafile'},...
    {'Other channel groups are'},...
    {'''EEG1010''    with approximately 90 electrodes'},...
    {'''EEG1005''    with approximately 350 electrodes'},...
    {'''EEGCHWILLA'' for Dorothee Chwilla''s electrode caps (used at the NICI)'},...
    {'''EEGBHAM''    for the 128 channel EEG system used in Birmingham'},...
    {'''EEGREF''     for mastoid and ear electrodes (M1, M2, LM, RM, A1, A2)'},...
    {'''MZ''         for MEG central'},...
    {'''ML''         for MEG left'},...
    {'''MR''         for MEG right'},...
    {'''MLx'', ''MRx'' and ''MZx'' with x=C,F,O,P,T for left/right central, frontal, occipital, parietal and temporal'},...
    {'You can also exclude channels or channel groups using the following syntax'},...
    {'{''all'', ''-POz'', ''-Fp1'', ''-EOG''}'}];
% ------------ Filters ------------


padding = struct('type','entry','name','Filter padding','tag','padding',...
    'strtype', 'n', 'val', {{0}}, 'num', [1 1], 'help', {{'Length to which the trials are padded for filtering'}});

% Low pass

lpfreq  = struct('type','entry','name','LP frequency','tag','lpfreq',...
    'strtype', 'r', 'num', [1 1], 'help', {{'Lowpass  frequency in Hz'}});

lpfiltord  = struct('type','entry','name','LP order','tag','lpfiltord',...
    'strtype', 'n', 'val', {{6}}, 'num', [1 1], 'help', {{'Lowpass  filter order'}});

lpfilttype = struct('type','menu','name','LP filter type','tag','lpfilttype',...
    'labels', {{'Butterworth', 'FIR'}}, 'values', {{'but', 'fir'}}, 'val', {{'but'}},  'help', {{'Type of the lowpass filter'}});

lpfiltdir = struct('type','menu','name','LP filter direction','tag','lpfiltdir',...
    'labels', {{'zero phase', 'forward', 'backward'}}, 'values', {{'twopass', 'onepass', 'onepass-reverse'}},...
    'val', {{'twopass'}},  'help', {{'Direction of the lowpass filter'}});

lpfilter = struct('type','branch','name','Lowpass filter','tag','lpfilter',...
    'val', {{lpfreq,lpfiltord, lpfilttype, lpfiltdir}} , 'help', {{'Should a lowpass filter be used'}});

% High pass

hpfreq  = struct('type','entry','name','HP frequency','tag','hpfreq',...
    'strtype', 'r', 'num', [1 1], 'help', {{'Highpass  frequency in Hz'}});

hpfiltord  = struct('type','entry','name','HP order','tag','hpfiltord',...
    'strtype', 'n', 'val', {{6}}, 'num', [1 1], 'help', {{'Highpass  filter order'}});

hpfilttype = struct('type','menu','name','HP filter type','tag','hpfilttype',...
    'labels', {{'Butterworth', 'FIR'}}, 'values', {{'but', 'fir'}}, 'val', {{'but'}},  'help', {{'Type of the highpass filter'}});

hpfiltdir = struct('type','menu','name','HP filter direction','tag','hpfiltdir',...
    'labels', {{'zero phase', 'forward', 'backward'}}, 'values', {{'twopass', 'onepass', 'onepass-reverse'}},...
    'val', {{'twopass'}},  'help', {{'Direction of the highpass filter'}});

hpfilter = struct('type','branch','name','Highpass filter','tag','hpfilter',...
    'val', {{hpfreq,hpfiltord, hpfilttype, hpfiltdir}} , 'help', {{'Should a highpass filter be used'}});


% Band pass

bpfreq  = struct('type','entry','name','BP frequency','tag','bpfreq',...
    'strtype', 'r', 'num', [1 2], 'help', {{'Bandpass  frequency in Hz'}});

bpfiltord  = struct('type','entry','name','BP order','tag','bpfiltord',...
    'strtype', 'n', 'val', {{4}}, 'num', [1 1], 'help', {{'Bandpass filter order'}});

bpfilttype = struct('type','menu','name','BP filter type','tag','bpfilttype',...
    'labels', {{'Butterworth', 'FIR'}}, 'values', {{'but', 'fir'}}, 'val', {{'but'}},  'help', {{'Type of the bandpass filter'}});

bpfiltdir = struct('type','menu','name','BP filter direction','tag','bpfiltdir',...
    'labels', {{'zero phase', 'forward', 'backward'}}, 'values', {{'twopass', 'onepass', 'onepass-reverse'}},...
    'val', {{'twopass'}},  'help', {{'Direction of the bandpass filter'}});

bpfilter = struct('type','branch','name','Bandpass filter','tag','bpfilter',...
    'val', {{bpfreq,bpfiltord, bpfilttype, bpfiltdir}} , 'help', {{'Should a bandpass filter be used'}});

%  Notch


lnfreq  = struct('type','entry','name','Notch frequency','tag','lnfreq',...
    'strtype', 'r', 'val', {{50}}, 'num', [1 1], 'help', {{'Notch  frequency in Hz'}});

lnfiltord  = struct('type','entry','name','Notch filter order','tag','lnfiltord',...
    'strtype', 'n', 'val', {{4}}, 'num', [1 1], 'help', {{'Notch  filter order'}});

lnfilter = struct('type','branch','name','Line noise filter','tag','lnfilter',...
    'val', {{lnfreq, lnfiltord}} , 'help', {{'Should a notch filter be used'}});

% DFT

dftfreq  = struct('type','entry','name','DFT frequencies','tag','dftfreq',...
    'strtype', 'r', 'val', {{[50 100 150]}}, 'num', [1 3], 'help', {{'Line noise and harmonics  frequencies in Hz'}});

dftfilter = struct('type','branch','name','DFT filter','tag','dftfilter',...
    'val', {{dftfreq}} , 'help', {{'Should a DFT filter be used'}});

% Median

medianfiltord  = struct('type','entry','name','Median filter order','tag','medianfiltord',...
    'strtype', 'n', 'val', {{9}}, 'num', [1 1], 'help', {{'Median filter order'}});

medianfilter = struct('type','branch','name','Median Filter','tag','medianfilter', 'val', {{medianfiltord}});

% Detrending and baseline correction

detrend = struct('type','const','name','Detrend','tag','detrend', 'val', {{'yes'}},  'help', {{'Detrend the data'}});

blcwindow  = struct('type','entry','name','Baseline window','tag','blcwindow',...
    'strtype', 'r', 'num', [1 2], 'help', {{'[begin end] in seconds, the default is the complete trial'}});

blc = struct('type','branch','name','Baseline correction','tag','blc', 'val', {{blcwindow}});

% Other
hilbert = struct('type','const','name','Hilbert transform','tag','hilbert', 'val', {{'yes'}},  'help', {{'Use Hilbert transform'}});

rectify  = struct('type','const','name','Rectify data','tag','rectify', 'val', {{'yes'}},  'help', {{'Rectify the data'}});

filters = struct('type','repeat','name','Filters and transformations','tag','filters', 'values',...
    {{lpfilter, hpfilter, bpfilter, lnfilter, dftfilter, medianfilter, detrend, blc, hilbert, rectify}}, 'check', @check_filter);


chansubset = struct('type', 'branch', 'name', 'Channel subset', 'tag', 'chansubset', 'val', {{channel, filters}});

channelsets = struct('type', 'repeat', 'name', 'Channels', 'tag', 'channelsets', 'values', {{chansubset}});

precision = struct('type','menu','name','Number precision','tag','precision',...
    'labels', {{'Single', 'Double'}}, 'values', {{'single', 'double'}}, 'val', {{'double'}},  'help', {{'Number precision for data processing'}});


% Save settings

ctf = struct('type','files','name','Name of the channel template file','tag','ctf',...
    'filter','mat','num',1, 'val', {{'empty.mat'}}, 'dir', fullfile(spm('dir'), 'EEGtemplates'), ...
    'help',{{'Select the channel template file (if available)'}});

savename  = struct('type','entry','name','Name of the output file','tag','savename',...
    'strtype', 's', 'num', [1 1], 'help', {{'Name of the file in which the data will be saved.'}});

sub_savesettings = struct('type','branch','name','Save settings','tag','sub_savesettings',...
    'val',{{savename, ctf}});

S   = struct('type','branch','name','MEEG Preprocessing','tag','meeg',...
    'val',{{dataset, trials, forcecont, channelsets, precision, sub_savesettings}},'prog',@meeg_preprocess,'modality',{{'EEG'}},...
    'help',{{'Preprocessing of MEG/EEG data'}});

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION runs the Fieldtrip preprocessing based on the config tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meeg_preprocess(S)

addpath(fullfile(spm('dir'), 'MEEGpreprocessing'), ...
    fullfile(spm('dir'), 'MEEGfileio'));

cfg.dataset=cell2mat(S.dataset);

switch cell2mat(fieldnames(S.trials))
    case 'sub_trialdef'
        cfg.trialdef=jobtree2cfg(S.trials.sub_trialdef);
        % In case the event value is numeric tries to convert it into a number
        try cfg.trialdef.eventvalue=str2num(cfg.trialdef.eventvalue); catch end
        cfg.trl=trialfun_spm(cfg);
    case 'trlfile'
        trl=load(cell2mat(S.trials.trlfile));
        if isstruct(trl)
            trl=trl.trl;
        end
        if size(trl, 2)<3
            error('TRL matrix should have at least 3 columns');
        elseif size(trl, 2)==3
            warning('For SPM the fourth column of TRL matrix should contain event codes. Adding 1s');
            trl=[trl ones(size(trl,1), 1)];
        end
        cfg.trl=trl;
    case 'loadcont'
        cfg.trialdef.triallength = Inf;
        cfg.trialdef.ntrials = 1;
        cfg.trl=trialfun_spm(cfg);
end

cfg.precision=S.precision;


cfg.forcecont=strcmpi(S.forcecont, 'yes');

% Runs separate preprocessing for different channel groups
data={};
for ind= 1:length(S.chansubset)
    if iscell(S.chansubset(ind))
        chansettings=S.chansubset{ind};
    else
        chansettings=S.chansubset(ind);
    end

    cfg1=cfg;
    cfg1.channel=chansettings.channel;

    for ind2=1:length(chansettings.filters)
        cfg1=mergestruct(cfg1, jobtree2cfg(chansettings.filters{ind2}));

        % This is to add 'yes' for the filter settings that appear in the
        % cfg tree
        cfg1=setfield(cfg1, cell2mat(fieldnames(chansettings.filters{ind2})), 'yes');
    end

    data{ind} = spm_eeg_preprocessing(cfg1); %spm_eeg_preprocessing(cfg1);
end

if length(data)>1
    data=ft_appenddata([], data{:});
else
    data=data{1};
end

[junk, ctfname] = fileparts(S.sub_savesettings.ctf{:});

ft2spm(data, S.sub_savesettings.savename, [ctfname '.mat']);


%%

function cfg=jobtree2cfg(jobtree)
% Converts the GUI tree to FT cfg struct
fields_list=fieldnames(jobtree);
cfg=[];
for ind =1:length(fields_list)
    if isstruct(getfield(jobtree, fields_list{ind}))
        % If there is a substruct in the tree calls itself recursively
        subcfg=jobtree2cfg(getfield(jobtree, fields_list{ind}));
        if  strncmp(fields_list{ind}, 'sub_', 4)
            % This makes it possible to keep substructs in the cfg
            % if their names are preceded by sub_ prefix
            cfg=setfield(cfg, fields_list{ind}(5:end), subcfg);
        else
            cfg=mergestruct(cfg, subcfg);
        end
    else
        value=getfield(jobtree, fields_list{ind});
        if iscell(value)
            if length(value)>1
                for ind2=1:length(value)
                    cfg=mergestruct(cfg, jobtree2cfg(value{ind2}));
                end
            else
                if ~isempty(value)
                    value=value{:};
                end
                cfg=setfield(cfg, fields_list{ind}, value);
            end
        else
            cfg=setfield(cfg, fields_list{ind}, value);
        end
    end
end


function str = check_filter(job)
fields={};
for n=1:length(job)
    fields=[fields fieldnames(job{n})];
end
if length(fields)~=length(unique(fields));
    str = 'Only one instance of each transformation is allowed.';
else
    str = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that merges the fields of two structs into one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function smerged=mergestruct(struct1, struct2)

if isempty(struct1)
    smerged=struct2;
    return;
end
if isempty(struct2)
    smerged=struct1;
    return;
end

names1 = fieldnames(struct1);
names2 = fieldnames(struct2);

if ~isempty(intersect(names1, names2))
    error('Structures with common fields cannot be merged');
end

data1=struct2cell(struct1);
data2=struct2cell(struct2);

smerged=cell2struct([data1; data2], [names1; names2], 1);
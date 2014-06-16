function spm_eeg_var_measures
% Function for computing Fourier coherence using Fieldtrip and VAR based directed measures
% using SPM's spectral toolbox, developed by Will Penny.
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL

% Vladimir Litvak
% $Id: spm_eeg_var_measures.m 6047 2014-06-16 11:11:22Z vladimir $

[Finter,Fgraph] = spm('FnUIsetup','MEEGtoools VAR measures', 0);

if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox','spectral'));
    addpath(fullfile(spm('Dir'),'toolbox','mixture'));
end

D = spm_eeg_load;

chan = spm_input('What channels?','+1','LFP|EEG|GUI');

if strcmp('GUI', chan)
    [selection, ok]= listdlg('ListString', D.chanlabels, 'SelectionMode', 'multiple' ,'Name', 'Select channels' , 'ListSize', [400 300]);
    if ~ok
        return;
    end
else
   selection = strmatch(chan, D.chantype);
end

if length(selection)>1
    chan = D.chanlabels(selection);
else
    error('Not enough channels');
end

if strcmp(D.type, 'continuous') 
    trldur = spm_input('Input trial length in sec)', '+1', 'r', '2', 1);
    
    cfg = [];
    cfg.trackcallinfo  = 'no';
    cfg.dataset = fullfile(D.path, D.fname);
    cfg.channel= chan;
    cfg.trl = 1:round(trldur*D.fsample):D.nsamples;
    cfg.trl = [cfg.trl(1:(end-1))' cfg.trl(2:end)'-1 zeros(length(cfg.trl)-1, 1)];
    data = ft_preprocessing(cfg);   
elseif D.ntrials < 3 || (D.nsamples/D.fsample) < 1
    error('Expecting continuous data or epoched files with multiple trials at least 1 sec long');
else % ============ Select the data and convert to Fieldtrip struct    
    clb = D.condlist;
    
    if numel(clb) > 1
        
        [selection, ok]= listdlg('ListString', clb, 'SelectionMode', 'multiple' ,'Name', 'Select conditions' , 'ListSize', [400 300]);
        
        if ~ok
            return;
        end
    else
        selection = 1;
    end
    
    ind = D.indtrial(clb(selection), 'GOOD');
    
    
    if isempty(ind)
        error('No valid trials found');
    end
    
    data = D.ftraw;
    data.trial = data.trial(ind);
    data.time =  data.time(ind);
end
%%

cfg=[];
cfg.channel = chan;
cfg.resamplefs = 250;
cfg.detrend = 'yes';
cfg.trackcallinfo  = 'no';
data = ft_resampledata(cfg, data);
%%
cfg =[];
cfg.channel = chan;
cfg.keeptrials = 'yes';
cfg.trackcallinfo  = 'no';
data= ft_timelockanalysis(cfg, data);

Ntrials=size(data.trial,1);
Ntime=length(data.time);

%%
cfg = [];
cfg.output ='powandcsd';
cfg.keeptrials = 'yes';
cfg.keeptapers='no';
cfg.taper = 'dpss';
cfg.method = 'mtmfft';
cfg.foilim     = [0 100]; % Frequency range
cfg.tapsmofrq = 1; % Frequency resolution
cfg.trackcallinfo  = 'no';
%
inp = ft_freqanalysis(cfg, data);

cfg1 = [];
cfg1.method = 'coh';
coh = ft_connectivityanalysis(cfg1, inp);
%%
% Defines how trials are shifted for shift-predictors
shift=[2:Ntrials 1];

% Compute the shift predictor for coherence
scoh=coh;
scoh.cohspctrm=zeros(size(scoh.cohspctrm));
for c=1:length(data.label)
    sdata=data;
    sdata.trial(:,c,:)=data.trial(shift,c,:);
    cfg.channelcmb = {data.label{c}, 'all'};
    inp = ft_freqanalysis(cfg, sdata);
    sscoh = ft_connectivityanalysis(cfg1, inp);
    for i=1:size(sscoh.labelcmb, 1)
        ind=[intersect(strmatch(sscoh.labelcmb(i,1),scoh.labelcmb(:,1),'exact'), ...
            strmatch(sscoh.labelcmb(i,2),scoh.labelcmb(:,2),'exact'))...
            intersect(strmatch(sscoh.labelcmb(i,1),scoh.labelcmb(:,2),'exact'), ...
            strmatch(sscoh.labelcmb(i,2),scoh.labelcmb(:,1),'exact'))];
        scoh.cohspctrm(ind, :)=sscoh.cohspctrm(i, :);
    end
end

% Plot power and coherence

figure('Name', 'Fourier power and coherence');
for i=1:1:length(data.label)
    for j=1:1:length(data.label)
        subplot(length(data.label),length(data.label), sub2ind([length(data.label) length(data.label)], i, j));
        if (i~=j)            
            ind=[intersect(strmatch(data.label{i},coh.labelcmb(:,1),'exact'), ...
                strmatch(data.label{j},coh.labelcmb(:,2),'exact'))...
                intersect(strmatch(data.label{i},coh.labelcmb(:,2),'exact'), ...
                strmatch(data.label{j},coh.labelcmb(:,1),'exact'))];
            plot(coh.freq, squeeze(coh.cohspctrm(ind,:)), 'b');
            hold on
            plot(coh.freq, squeeze(scoh.cohspctrm(ind,:)), 'r');
            ylim([0 1]);
            xlim([0 100]);
            title([data.label{i} '->' data.label{j}]);
        else
            plot(inp.freq, log(squeeze(mean(inp.powspctrm(:, i, :),1))), 'b');
            xlim([0 100]);
            title(data.label{i});
        end
    end
end

%% This is the code for determining optimal model order

if spm_input('Check model order?','+1','yes|no',[1 0]);

    F=[];
    testdata=squeeze(data.trial(1,:,:))';
    max_p=15; % Maximal order to check
    for p=1:max_p,
        mar=spm_mar(testdata,p);
        F(p)=mar.fm;
        disp(sprintf('Model order p=%d, Evidence = %1.2f',p,F(p)));
    end
    Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
    plot([1:max_p],F);
    title('all regions');
    xlabel('Model order,p');
    ylabel('Log Evidence');
end



%%
% Compute the MAR model 
ns=1/mean(diff(data.time)); % Sample rate
p= spm_input('Input model order', '+1', 'r', '', 1); % Order of MAR model
freqs=[1:100]; % Frequencies to evaluate spectral quantities at

dircoh=[];
sdircoh=[];

% Here the specific directional measure can be selected
measure= spm_input('Select VAR measure?','+1','pdc|pve|dtf|C'); 

for tr=1:Ntrials

    ctrial=prestd(squeeze(data.trial(tr,:,:)));

    mar = spm_mar(ctrial,p);
    mar = spm_mar_spectra(mar,freqs,ns);
    dircoh(tr, :, :, :)=squeeze(getfield(mar, measure));

    strial=prestd(squeeze(data.trial(shift(tr),:,:)));

    % Compute the shift predictors for each channel separately
    for c=1:length(data.label)

        combtrial=ctrial;
        combtrial(:, c)=strial(:, c);
        mar = spm_mar (combtrial ,p);
        mar = spm_mar_spectra (mar,freqs,ns);
        sdircoh(c, tr, :, :, :)=squeeze(getfield(mar, measure));
    end
    disp(sprintf('Trial %d out of %d',tr,Ntrials));
end

%%  Put the results back in FieldTrip data structures

dummy=[]; % The data
dummy.freq=freqs;
dummy.dimord='rpt_chan_freq';
dummy.trial=[];
dummy.label={};

sdummy=dummy; % The shift predictor

for i=1:length(data.label)
    for j=1:length(data.label)
        if(i~=j)
            dummy.trial=cat(2, dummy.trial, permute(shiftdim((dircoh(:, :, j, i)), -1), [2 1 3]));
            sdummy.trial=cat(2, sdummy.trial, permute(shiftdim(squeeze(sdircoh(j,:, :, j, i)), -1), [2 1 3]));
            dummy.label=[dummy.label; {[data.label{i} '->' data.label{j}]}];
            sdummy.label=[sdummy.label; {[data.label{i} '->' data.label{j}]}];
        end
    end
end

%% Run nonparametric comparison between the data and the shift predictor

Ntrials1=size(dummy.trial,1);
Ntrials2=size(sdummy.trial,1);

cfg=[];
cfg.frequency   = [1 100];
cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.tail = 1;
cfg.alpha = 0.05;
cfg.numrandomization = 1000;

cfg.correctm = spm_input('Choose MC correction', 1, 'm', 'no|max|cluster|bonferoni|holms|fdr');

cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail=1;


design = zeros(1,Ntrials1+Ntrials2);
design(1,1:Ntrials1) = 1;
design(1,(Ntrials1+1):end)=2;


cfg.design   = design;   % design matrix
cfg.ivar  = 1;  % number or list with indices, independent variable(s)

cfg.neighbours=[];
for c=1:length(dummy.label)
    cfg.neighbours(c).label = dummy.label{c};
    cfg.neighbours(c).neighblabel ={};
end

cfg.parameter='trial';
%%
stats= ft_freqstatistics(cfg, dummy, sdummy);

%% Plot the results of the MAR analysis

mdircoh=squeeze(mean(dummy.trial,1));
smdircoh=squeeze(mean(sdummy.trial,1));
figure('Name', ['VAR-' upper(measure)]);
for i=1:1:length(data.label)
    for j=1:1:length(data.label)
        if (i~=j)
            subplot(length(data.label),length(data.label), sub2ind([length(data.label) length(data.label)], i, j));
            maskind=spm_match_str(dummy.label, [data.label{i} '->' data.label{j}]);
            plot(freqs, squeeze(mdircoh(maskind,:)));
            sigind=find(squeeze(stats.mask(maskind,:)));
            hold on
            plot(freqs(sigind), squeeze(mdircoh(maskind, sigind)), '*');
            plot(freqs, squeeze(smdircoh(maskind,:)), 'r');
            ylim([0 max(max(max(mdircoh)), max(max(smdircoh)))]);
            xlim([0 100]);
            title([data.label{i} '->' data.label{j}]);
        end
    end
end

%%
function X=prestd(X)

if size(X, 2)>size(X,1)
    X=X';
end

X=X-repmat(mean(X), size(X,1),1);
X=X./repmat(std(X), size(X,1),1);

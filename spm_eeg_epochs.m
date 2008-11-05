function D = spm_eeg_epochs(S)
% function used for epoching continuous EEG/MEG data
% FORMAT D = spm_eeg_epochs(S)
%
% S  - filename or input struct (optional)
% (optional) fields of S:
% S.D         - filename of EEG mat-file with continuous data
% 
% Either (to use a ready-made trial definition): 
% S.epochinfo.trl - Nx2 or Nx3 matrix (N - number of trials) [start end offset]
% S.epochinfo.conditionlabels - one label or cell array of N labels 
% S.epochinfo.padding - in sec - the additional time period around each trial
%               for which the events are saved with the trial (to let the
%               user keep and use for analysis events which are outside
% 
% or (to define trials using (spm_eeg_definetrial)
%
% S.pretrig - pre-trigger time in ms
% S.posttrig - post-trigger time in ms. 
% S.trialdef - structure array for trial definition with fields
%       S.trialdef.conditionlabel - string label for the condition
%       S.trialdef.eventtype  - string
%       S.trialdef.eventvalue  - string, numeric or empty
% S.review - 1 - review individual trials after selection
% S.save   - 1 - save trial definition
%
% Output:
% D         - EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_epochs extracts single trials from continuous EEG/MEG data. The
% length of an epoch is determined by the samples before and after stimulus
% presentation. One can limit the extracted trials to specific trial types.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_epochs.m 2446 2008-11-05 16:05:14Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG epoching setup',0);

if nargin == 0
    S =[];
end

try
    D = S.D;
catch
    D = spm_select(1, 'mat', 'Select M/EEG mat file');
end

try
    D = spm_eeg_load(D);
catch
    error(sprintf('Trouble reading file %s', D));
end

if ~strncmpi(D.type, 'cont', 4)
    warning('The file must contain continuous data.');
    return
end


% first case: deftrials (default for GUI call)
if isfield(S, 'trialdef') || nargin == 0
    S_definetrial = [];
    
    if isfield(S, 'pretrig')
        S_definetrial.pretrig = S.pretrig;
    end
      
    if isfield(S, 'posttrig')
        S_definetrial.posttrig = S.posttrig;
    end
    
    if isfield(S, 'trialdef')
        S_definetrial.trialdef = S.trialdef;
    end
    
    if isfield(S, 'review')
        S_definetrial.review = S.review;
    end
    
    if isfield(S, 'save')
        S_definetrial.save = S.save;
    end
    
    S_definetrial.event = D.events;
    
    S_definetrial.fsample = D.fsample;
    
    S_definetrial.timeonset = D.timeonset;
    
    [epochinfo.trl, epochinfo.conditionlabels, S] = spm_eeg_definetrial(S_definetrial);

% second case: epochinfo (trlfile and trl)
else
    try
        epochinfo.trl = S.epochinfo.trl;
        epochinfo.conditionlabels = S.epochinfo.conditionlabels;
    catch
        try
            epochinfo.trlfile = S.epochinfo.trlfile;
        catch
            epochinfo.trlfile = spm_select(1, '\.mat$', 'Select a trial definition file');
        end
        try
            epochinfo.trl = getfield(load(S.epochinfo.trlfile, 'trl'), 'trl');
            epochinfo.conditionlabels = getfield(load(epochinfo.trlfile, 'conditionlabels'), 'conditionlabels');
        catch
            error('Trouble reading trl file.');
        end
    end
   
end

trl = epochinfo.trl;
conditionlabels = epochinfo.conditionlabels;
S.D = fullfile(D.path, D.fname); % history

try
    epochinfo.padding = S.epochinfo.padding;
catch
    epochinfo.padding = 0;
    % for history
    S.epochinfo.padding = epochinfo.padding;
end


% checks on input
if size(trl, 2) >= 3
    timeOnset = unique(trl(:, 3))./D.fsample;
    trl = trl(:, 1:2);
else
    timeOnset = 0;
end

if length(timeOnset) > 1
    error('All trials should have identical baseline');
end

nsampl = unique(round(diff(trl, [], 2)))+1;
if length(nsampl) > 1 || nsampl<1
    error('All trials should have identical and positive lengths');
end

spm('Pointer', 'Watch'); drawnow;

inbounds = (trl(:,1)>1 & trl(:, 2)<=D.nsamples);

rejected = find(~inbounds);
rejected = rejected(:)';

if ~isempty(rejected)
    trl = trl(find(inbounds), :);
    warning([D.fname ': Events ' num2str(rejected) ' not extracted - out of bounds']);
end

ntrial = size(trl, 1);

% generate new meeg object with new filenames
Dnew = clone(D, ['e' fnamedat(D)], [D.nchannels nsampl, ntrial]);

spm_progress_bar('Init', ntrial, 'Events read'); drawnow;
if ntrial > 100, Ibar = floor(linspace(1, ntrial, 100));
else Ibar = [1:ntrial]; end

for i = 1:ntrial

    d = D(:, trl(i, 1):trl(i, 2), 1);

    Dnew(:, :, i) = d;

    Dnew = events(Dnew, i, select_events(D.events, ...
        [trl(i, 1)/D.fsample-epochinfo.padding  trl(i, 2)/D.fsample+epochinfo.padding]));

    if ismember(i, Ibar)
        spm_progress_bar('Set', i);
        drawnow;
    end
end

Dnew = conditions(Dnew, [], conditionlabels);

% The conditions will later be sorted in the original order they were defined.
if isfield(S, 'trialdef')
    Dnew = condlist(Dnew, {S.trialdef(:).conditionlabel});
end

Dnew = trialonset(Dnew, [], trl(:, 1)./D.fsample+D.trialonset);
Dnew = timeonset(Dnew, timeOnset);
Dnew = type(Dnew, 'single');



Dnew = spm_eeg_bc(Dnew, [time(Dnew, 1, 'ms') 0]);

% history (remove first some redundant stuff potentially put in by
% spm_eeg_definetrial)
if isfield(S, 'event')
    S = rmfield(S, 'event');
end

D = Dnew;
D = D.history('spm_eeg_epochs', S);
save(D);

spm_progress_bar('Clear');

spm('Pointer', 'Arrow');

function event = select_events(event, timeseg)
% Utility function to select events according to time segment
% FORMAT event = select_events(event, timeseg)

if ~isempty(event)
    [time ind] = sort([event(:).time]);

    selectind = ind(time>=timeseg(1) & time<=timeseg(2));

    event = event(selectind);
end

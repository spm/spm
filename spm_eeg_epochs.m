function D = spm_eeg_epochs(S)
% function used for epoching continuous EEG/MEG data
% FORMAT D = spm_eeg_epochs(S)
%
% S  - filename or input struct (optional)
% (optional) fields of S:
% S.D         - filename of EEG mat-file with continuous data
% S.trl - Nx2 or Nx3 matrix (N - number of trials) [start end offset]
% S.conditionlabels - one label or cell array of N labels 
% S.eventpadding - in sec - the additional time period around each trial
%               for which the events are saved with the trial (to let the
%               user keep and use for analysis events which are outside
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
% $Id: spm_eeg_epochs.m 1254 2008-03-27 18:41:42Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG epoching setup',0);

if nargin == 0
    S =[];
end

if ischar(S)
    temp = S;
    S=[];
    S.D = temp;
end

if ~isfield(S, 'eventpadding'),    S.eventpadding = 0;    end

if isfield(S, 'D')
    D = S.D;
else
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
end


if isa(D, 'char')
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

if ntrials(D)>1
    warning('The file is already epoched');
    return;
end

if ~(isfield(S, 'trl') & isfield(S, 'conditionlabels'))
    S.event = D.events;
    if isempty(S.event)
        S.trlfile = spm_select(1, '\.mat$', 'Select a trial definition file');
        S.trl = getfield(load(S.trlfile, 'trl'), 'trl');
        S.conditionlabels = getfield(load(S.trlfile, 'conditionlabels'), 'conditionlabels');
    end
    S.fsample = D.fsample;
    S.timeonset = D.trialonset;
    [S.trl, S.conditionlabels] = spm_eeg_definetrial(S);
end

trl = S.trl;
conditionlabels = S.conditionlabels;

if size(trl, 2) >= 3
    timeOnset = unique(trl(:, 3))./D.fsample;
    trl = trl(:, 1:2);
else
    timeOnset = 0;
end

if length(timeOnset) > 1
    error('All trials should have identical baseline');
end

nsampl = unique(diff(trl, [], 2))+1;
if length(nsampl) > 1 || nsampl<1
    error('All trials should have identical and positive lengths');
end

spm('Pointer', 'Watch'); drawnow;

inbounds = (trl(:,1)>1 & trl(:, 2)<=D.nsamples);

rejected = find(~inbounds);

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
        [trl(i, 1)/D.fsample-S.eventpadding  trl(i, 2)/D.fsample+S.eventpadding]));

    if ismember(i, Ibar)
        spm_progress_bar('Set', i);
        drawnow;
    end
end

Dnew = conditions(Dnew, [], conditionlabels);
Dnew = trialonset(Dnew, [], trl(i, 1)./D.fsample+D.trialonset);
Dnew = timeonset(Dnew, timeOnset);

save(Dnew);

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

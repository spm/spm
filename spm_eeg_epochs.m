function D = spm_eeg_epochs(S)
% Epoching continuous M/EEG data
% FORMAT D = spm_eeg_epochs(S)
%
% S                     - input structure (optional)
% (optional) fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file with
%                         continuous data
%   S.bc                - baseline-correct the data (1 - yes, 0 - no).
%
% Either (to use a ready-made trial definition):
%   S.epochinfo.trl             - Nx2 or Nx3 matrix (N - number of trials)
%                                 [start end offset]
%   S.epochinfo.conditionlabels - one label or cell array of N labels
%   S.epochinfo.padding         - the additional time period around each
%                                 trial for which the events are saved with
%                                 the trial (to let the user keep and use
%                                 for analysis events which are outside) [in ms]
%
% Or (to define trials using (spm_eeg_definetrial)):
%   S.pretrig           - pre-trigger time [in ms]
%   S.posttrig          - post-trigger time [in ms]
%   S.trialdef          - structure array for trial definition with fields
%     S.trialdef.conditionlabel - string label for the condition
%     S.trialdef.eventtype      - string
%     S.trialdef.eventvalue     - string, numeric or empty
%
%   S.reviewtrials      - review individual trials after selection
%   S.save              - save trial definition
%
% Output:
% D                     - MEEG object (also written on disk)
%__________________________________________________________________________
%
% spm_eeg_epochs extracts single trials from continuous EEG/MEG data. The
% length of an epoch is determined by the samples before and after stimulus
% presentation. One can limit the extracted trials to specific trial types.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_epochs.m 4244 2011-03-14 13:33:01Z vladimir $

SVNrev = '$Rev: 4244 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG epoching'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D   = spm_eeg_load(D);
S.D = fullfile(D.path, D.fname);

%-Check that the input file contains continuous data
%--------------------------------------------------------------------------
if ~strncmpi(D.type, 'cont', 4)
    error('The file must contain continuous data.');
end


if ~isfield(S, 'bc')
    S.bc = 1;
    % S.bc = spm_input('Subtract baseline?','+1','yes|no',[1 0], 1);
end

%-First case: deftrials (default for GUI call)
%--------------------------------------------------------------------------
if isfield(S, 'trialdef') || nargin == 0

    if isfield(S, 'pretrig')
        S_definetrial.pretrig = S.pretrig;
    end

    if isfield(S, 'posttrig')
        S_definetrial.posttrig = S.posttrig;
    end

    if isfield(S, 'trialdef')
        S_definetrial.trialdef = S.trialdef;
    end

    if isfield(S, 'reviewtrials')
        S_definetrial.reviewtrials = S.reviewtrials;
    end

    if isfield(S, 'save')
        S_definetrial.save = S.save;
    end

    S_definetrial.D     = S.D;
    
    S_definetrial.event = D.events;

    S_definetrial.fsample = D.fsample;

    S_definetrial.timeonset = D.timeonset;

    S_definetrial.bc = S.bc;

    [epochinfo.trl, epochinfo.conditionlabels, S] = spm_eeg_definetrial(S_definetrial);

    %-Second case: epochinfo (trlfile and trl)
    %--------------------------------------------------------------------------
else
    try
        epochinfo.trl = S.epochinfo.trl;
        epochinfo.conditionlabels = S.epochinfo.conditionlabels;
    catch
        try
            epochinfo.trlfile = S.epochinfo.trlfile;
        catch
            epochinfo.trlfile = spm_select(1, 'mat', 'Select a trial definition file');
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

if numel(conditionlabels) == 1
   conditionlabels = repmat(conditionlabels, 1, size(trl, 1));
end

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

inbounds = (trl(:,1)>=1 & trl(:, 2)<=D.nsamples);

rejected = find(~inbounds);
rejected = rejected(:)';

if ~isempty(rejected)
    trl = trl(find(inbounds), :);
    warning([D.fname ': Events ' num2str(rejected) ' not extracted - out of bounds']);
end

ntrial = size(trl, 1);

%-Generate new MEEG object with new filenames
%--------------------------------------------------------------------------
Dnew = clone(D, ['e' fnamedat(D)], [D.nchannels nsampl, ntrial]);

%-Epoch data
%--------------------------------------------------------------------------
spm_progress_bar('Init', ntrial, 'Events read');
if ntrial > 100, Ibar = floor(linspace(1, ntrial, 100));
else Ibar = [1:ntrial]; end

for i = 1:ntrial

    d = D(:, trl(i, 1):trl(i, 2), 1);

    Dnew(:, :, i) = d;

    Dnew = events(Dnew, i, select_events(D.events, ...
        [trl(i, 1)/D.fsample-epochinfo.padding  trl(i, 2)/D.fsample+epochinfo.padding]));

    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

Dnew = conditions(Dnew, [], conditionlabels);

% The conditions will later be sorted in the original order they were defined.
if isfield(S, 'trialdef')
    Dnew = condlist(Dnew, {S.trialdef(:).conditionlabel});
end

Dnew = trialonset(Dnew, [], trl(:, 1)./D.fsample+D.trialonset);
Dnew = timeonset(Dnew, timeOnset);
Dnew = type(Dnew, 'single');

%-Perform baseline correction if there are negative time points
%--------------------------------------------------------------------------
if S.bc
    if time(Dnew, 1) < 0
        S1               = [];
        S1.D             = Dnew;
        S1.time          = [time(Dnew, 1, 'ms') 0];
        S1.save          = false;
        S1.updatehistory = false;
        Dnew             = spm_eeg_bc(S1);
    else
        warning('There was no baseline specified. The data is not baseline-corrected');
    end
end

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = Dnew;
% Remove some redundant stuff potentially put in by spm_eeg_definetrial
if isfield(S, 'event'), S = rmfield(S, 'event'); end
D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG epoching: done'); spm('Pointer','Arrow');


%==========================================================================
function event = select_events(event, timeseg)
% Utility function to select events according to time segment

if ~isempty(event)
    [time ind] = sort([event(:).time]);

    selectind = ind(time >= timeseg(1) & time <= timeseg(2));

    event = event(selectind);
end

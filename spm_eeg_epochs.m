function D = spm_eeg_epochs(S)
% Epoching continuous M/EEG data
% FORMAT D = spm_eeg_epochs(S)
%
% S                   - input structure 
%  fields of S:
%   S.D               - MEEG object or filename of M/EEG mat-file with
%                         continuous data
%   S.bc              - baseline-correct the data (1 - yes, 0 - no).
%
% Either (to use a ready-made trial definition):
%     S.trl            - [N x 3] trl matrix or name of the trial definition file
%                      containing 'trl' variable with such a matrix
%     S.conditionlabels - labels for the trials in the data [default: 'Undefined']
%
%  or
%
%     S.timewin         -  time window in PST ms
%     S.trialdef       - structure array for trial definition with fields
%       S.trialdef.conditionlabel - string label for the condition
%       S.trialdef.eventtype      - string
%       S.trialdef.eventvalue     - string, numeric or empty
%
%
%    S.eventpadding  -  (optional) the additional time period around each
%                                     trial for which the events are saved with
%                                     the trial (to let the user keep and use
%                                     for analysis events which are outside) [in ms]
%
%    S.prefix     - prefix for the output file (default - 'e')
%
%
% Output:
% D                     - MEEG object (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_epochs.m 5074 2012-11-23 12:18:26Z vladimir $

SVNrev = '$Rev: 5074 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG epoching'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix = 'e';           end
if ~isfield(S, 'bc'),           S.bc = 1;                 end
if ~isfield(S, 'eventpadding'), S.eventpadding = 0;       end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);


%-Check that the input file contains continuous data
%--------------------------------------------------------------------------
if ~isequal(D.type, 'continuous')
    error('The file must contain continuous data.');
end


if all(isfield(S, {'trialdef', 'timewin'}))
    S1 = [];
    S1.D = D;
    S1.reviewtrials = 0;
    S1.save = 0;
    
    if ischar(S.trialdef)
        S1.trialdef = getfield(load(S.trialdef), 'trialdef');
    else
        S1.trialdef = S.trialdef;
    end
    
    S1.timewin = S.timewin;
    
    [trl, conditionlabels] = spm_eeg_definetrial(S1);

elseif isfield(S, 'trl')
    if ischar(S.trl)
        trlfile = load(S.trl);
        trl = getfield(trlfile, 'trl');
        
        if isfield(trlfile, 'conditionlabels')
            conditionlabels = getfield(trlfile, 'conditionlabels');
        else
            conditionlabels = 'Undefined';
        end
    else
        trl = S.trl;
        if isfield(S, 'conditionlabels')
            conditionlabels = S.conditionlabels;
        else
            conditionlabels = 'Undefined';
        end
    end
else
    error('Invalid trial definition');
end
   
if ischar(conditionlabels)
    conditionlabels = {conditionlabels};
end

if numel(conditionlabels) == 1
   conditionlabels = repmat(conditionlabels, 1, size(trl, 1));
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
    trl = trl(inbounds, :);
    conditionlabels = conditionlabels(inbounds);
    warning([D.fname ': Events ' num2str(rejected) ' not extracted - out of bounds']);
end

ntrial = size(trl, 1);

%-Generate new MEEG object with new filenames
%--------------------------------------------------------------------------
Dnew = clone(D, [S.prefix fname(D)], [D.nchannels nsampl, ntrial]);

%-Baseline correction
%--------------------------------------------------------------------------
if S.bc
    if time(Dnew, 1) < 0
        bc = Dnew.indsample(0);
        chanbc = D.indchantype('Filtered');
    else
       bc = 0;
       warning('There was no baseline specified. The data is not baseline-corrected');
    end
end

%-Epoch data
%--------------------------------------------------------------------------
spm_progress_bar('Init', ntrial, 'Events read');
if ntrial > 100, Ibar = floor(linspace(1, ntrial, 100));
else Ibar = [1:ntrial]; end

for i = 1:ntrial

    d = D(:, trl(i, 1):trl(i, 2), 1);
    
    if bc
        mbaseline = mean(d(chanbc, 1:bc), 2);
        d(chanbc, :) = d(chanbc, :) - repmat(mbaseline, 1, size(d, 2));
    end

    Dnew(:, :, i) = d;

    Dnew = events(Dnew, i, select_events(D.events, ...
        [trl(i, 1)/D.fsample-S.eventpadding  trl(i, 2)/D.fsample+S.eventpadding]));

    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

Dnew = conditions(Dnew, ':', conditionlabels);

% The conditions will later be sorted in the original order they were defined.
if isfield(S, 'trialdef')
    Dnew = condlist(Dnew, {S.trialdef(:).conditionlabel});
end

Dnew = trialonset(Dnew, ':', trl(:, 1)./D.fsample+D.trialonset);
Dnew = timeonset(Dnew, timeOnset);
Dnew = type(Dnew, 'single');

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

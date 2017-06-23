function D = spm_eeg_epochs(S)
% Epoching continuous M/EEG data
% FORMAT D = spm_eeg_epochs(S)
%
% S                - input structure 
%  fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file with
%                         continuous data
%   S.bc                - baseline-correct the data [1: yes, 0: no]
%
% Either (to use a ready-made trial definition):
%
%     S.trl             - [N x 3] trl matrix or name of the trial definition
%                         file containing 'trl' variable with such a matrix
%
%     S.conditionlabels - labels for the trials in the data
%                         [default: 'Undefined']
%
%  or
%
%     S.timewin         - time window in PST ms
%
%     S.trialdef        - structure array for trial definition with fields
%       S.trialdef.conditionlabel - string label for the condition
%       S.trialdef.eventtype      - string
%       S.trialdef.eventvalue     - string, numeric or empty
%
%  or
%       
%    S.trialength       - length of arbitray trials to split the data into
%                         (in ms). This is useful e.g. for spectral
%                         analysis of steady state data
%
%    S.conditionlabels  - labels for the trials in the data
%                         [default: 'Undefined']
%
%    S.eventpadding     - (optional) the additional time period around each
%                         trial for which the events are saved with
%                         the trial (to let the user keep and use
%                         for analysis events which are outside) {in s}
%                         [default: 0]
%
%    S.prefix           - prefix for the output file [default: 'e']
%
%
% Output:
% D                     - MEEG object (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_epochs.m 7125 2017-06-23 09:49:29Z guillaume $

SVNrev = '$Rev: 7125 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG epoching'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

isTF = strncmpi(D.transformtype,'TF',2);

if isTF && isfield(S, 'bc') && S.bc
    sw = warning('off','backtrace');
    warning('Automatic baseline correction is not done for TF data. Use TF rescaling.');
    warning(sw);
    S.bc = false;
end

%-Input parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'prefix'),       S.prefix = 'e';     end
if ~isfield(S, 'eventpadding'), S.eventpadding = 0; end
if ~isfield(S, 'bc'),           S.bc = ~isTF;       end

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
        trl     = trlfile.trl;
        
        if isfield(trlfile, 'conditionlabels')
            conditionlabels = trlfile.conditionlabels;
        else
            conditionlabels = 'Undefined';
        end
    else
        trl     = S.trl;
        if isfield(S, 'conditionlabels')
            conditionlabels = S.conditionlabels;
        else
            conditionlabels = 'Undefined';
        end
    end
elseif isfield(S, 'trialength')
    trl = 1:round(1e-3*S.trialength*D.fsample):D.nsamples;
    trl = [trl(1:(end-1))' trl(2:end)' 0*trl(2:end)'];
    
    if isfield(S, 'conditionlabels')
        conditionlabels = S.conditionlabels;
    else
        conditionlabels = 'Undefined';
    end
else
    error('Invalid trial definition.');
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
    error('All trials should have identical baseline.');
end

if isempty(trl)
    error('No trials found.');
end

nsampl = unique(round(diff(trl, [], 2)))+1;
if length(nsampl) > 1 || nsampl<1
    error('All trials should have identical and positive lengths.');
end

inbounds = (trl(:,1)>=1 & trl(:, 2)<=D.nsamples);

rejected = find(~inbounds);
rejected = rejected(:)';

if ~isempty(rejected)
    trl = trl(inbounds, :);
    conditionlabels = conditionlabels(inbounds);
    warning([D.fname ': Events ' num2str(rejected) ' not extracted - out of bounds.']);
end

ntrial = size(trl, 1);

%-Generate new MEEG object with new filenames
%--------------------------------------------------------------------------
if isTF
    Dnew = clone(D, [S.prefix fname(D)], [D.nchannels, D.nfrequencies, nsampl, ntrial]);
else
    Dnew = clone(D, [S.prefix fname(D)], [D.nchannels, nsampl, ntrial]);
end

Dnew = timeonset(Dnew, timeOnset);
Dnew = type(Dnew, 'single');

%-Baseline correction
%--------------------------------------------------------------------------
if S.bc
    if time(Dnew, 1) < 0
        bc = Dnew.indsample(0);
        chanbc = D.indchantype('Filtered');
    elseif isfield(S, 'trialength')
        bc = Dnew.nsamples;
        chanbc = D.indchantype('Filtered');
    else
       S.bc = false;
       warning('There was no baseline specified. The data is not baseline-corrected.');
    end
end

%-Epoch data
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Baseline correction',num2str(S.bc));           %-#
fprintf('%-40s: %30s\n','Number of trials',num2str(ntrial));            %-#
spm_progress_bar('Init', ntrial, 'Trials completed');
if ntrial > 100, Ibar = floor(linspace(1, ntrial, 100));
else Ibar = [1:ntrial]; end

for i = 1:ntrial
    if isTF
        d = D(:, :, trl(i, 1):trl(i, 2), 1);
        Dnew(:, :, :, i) = d;
    else
        d = D(:, trl(i, 1):trl(i, 2), 1);
        
        if S.bc
            mbaseline = mean(d(chanbc, 1:bc), 2);
            d(chanbc, :) = d(chanbc, :) - repmat(mbaseline, 1, size(d, 2));
        end
        
        Dnew(:, :, i) = d;
    end
    
    Dnew = events(Dnew, i, select_events(D.events, ...
        D.trialonset+[trl(i, 1)/D.fsample-S.eventpadding  trl(i, 2)/D.fsample+S.eventpadding]));
    
    if any(Ibar == i), spm_progress_bar('Set', i); end
end

Dnew = conditions(Dnew, ':', conditionlabels);

% The conditions will later be sorted in the original order they were defined.
if isfield(S, 'trialdef')
    Dnew = condlist(Dnew, {S.trialdef(:).conditionlabel});
end

Dnew = trialonset(Dnew, ':', trl(:, 1)./D.fsample+D.trialonset);


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
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG epoching: done'); spm('Pointer','Arrow');


%==========================================================================
function event = select_events(event, timeseg)
% Utility function to select events according to time segment

if ~isempty(event)
    [time,ind] = sort([event(:).time]);

    selectind  = ind(time >= timeseg(1) & time <= timeseg(2));

    event      = event(selectind);
end

function D = spm_eeg_epochs(S)
% function used for epoching continuous EEG/MEG data
% FORMAT D = spm_eeg_epochs(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file with continuous data
% events    - struct with various entries:
%    start     - pre-stimulus start of epoch [ms]
%    stop      - post-stimulus end of epoch [ms]
%    types     - events to extract (vector of event types)
%    Inewlist  - switch (0/1) to have new list of event codes
%    Ec        - list of new event codes
%
% Output:
% D         - EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_epochs extracts single trials from continuous EEG/MEG data. The
% length of an epoch is determined by the samples before and after stimulus
% presentation. One can limit the extracted trials to specific trial types.
% Also, it is possible to re-number trial types, see above.
% Note that epoching includes a baseline correction of each single trial,
% i.e. a subtraction of the average pre-stimulus average from all time
% points.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_epochs.m 1247 2008-03-26 11:08:36Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG epoching setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_load(D);
catch
    error(sprintf('Trouble reading file %s', D));
end

try
    events.start = S.events.start;
catch
    events.start =...
        spm_input('start of epoch [ms]', '+1', 'r', '', 1);
end

try
    events.stop = S.events.stop;
catch
    events.stop = ...
        spm_input('end of epoch [ms]', '+1', 'r', '', 1);
end

% convert to struct
sD = struct(D);
ev = cat(2, sD.trials.events(:).value);
values = unique(ev);

try
    events.types = S.events.types;
catch
    disp(sprintf('Nr: Available events (#occurences): '))
    for i = 1:length(values)
        Nevents(i) = sum(values(i) == ev);
        disp(sprintf('%d:   %d (%d)', i, values(i), Nevents(i)))

    end
    events.types = ...
        spm_input('Event nrs. to epoch)', '+1', 'n', '', inf, length(values));
end



% returns indices of chosen trial types
selected = find(ismember(ev, values(events.types)));

try
    Inewlist = S.events.Inewlist;
catch
    Inewlist = spm_input('Read new event list?', '+1', 'yes|no', [1 0]);
end

if Inewlist
    try
        Ec = S.events.Ec;
    catch
        Ec = spm_input('Input event vector', '+1', 'w', [], length(values));
    end
end

try
    events.newlabels = S.events.newlabels;
catch
    events.newlabels = ...
        spm_input('New labels', '+1', 's+', '', length(events.types));
end

events.labels(events.types) = events.newlabels;

spm('Pointer', 'Watch'); drawnow;

% transform ms to samples
events.startsample = ceil(-events.start*D.fsample/1000);
events.stopsample = ceil(events.stop*D.fsample/1000);

if events.startsample >= events.stopsample
    error('Start time must be less than stop time.');
end

% two passes

% 1st pass: Count the number of trials to be epoched to know the
% dimensions of the resulting data array
ind = [];
t = round(cat(1, sD.trials.events(:).time)*D.fsample); % sample time
for i = selected
    if  t(i) - events.startsample < 1 || t(i) + events.stopsample > D.nsamples
        % skip this trial
        warning(sprintf('%s: Event %d not extracted because not enough sample points', D.fname, i));
    else
        ind = [ind i];
    end
end

% 2nd pass: do the epoching

nsamples = events.stopsample + events.startsample+1;
ntrials = length(ind);

% generate new meeg object with new filenames
Dnew = newdata(D, ['e' fnamedat(D)], [D.nchannels nsamples ntrials], D.dtype);

spm_progress_bar('Init', ntrials, 'Events read'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials, 100));
else Ibar = [1:ntrials]; end

for i = 1:ntrials

    k = ind(i);
    if Inewlist
        ev(i) = Ec(k);
    end

    d = D(:, t(k) - events.startsample : t(k) + events.stopsample, 1);

    % baseline subtraction
    d = d - repmat(mean(d(:, [1:abs(events.startsample)+1]), 2), 1, nsamples);

    Dnew(1:Dnew.nchannels, 1:nsamples, i) = d;

    if ismember(i, Ibar)
        spm_progress_bar('Set', i);
        drawnow;
    end
end

sDnew = struct(Dnew);
sDnew.Nsamples = nsamples;
trials = sDnew.trials;
sDnew = rmfield(sDnew, 'trials');

for i = 1: ntrials
    sDnew.trials(i) = struct('onset', t(ind(i))/D.fsample, 'label', num2str(events.labels{find(ev(ind(i)) == values)}));
end
Dnew = meeg(sDnew);

if Inewlist && D.ntrials ~= length(Ec)
    warning('Not all events in event list used!')
end

save(Dnew);

spm_progress_bar('Clear');

spm('Pointer', 'Arrow');

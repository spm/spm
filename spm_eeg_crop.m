function D = spm_eeg_crop(S)
% Reduce the data size by cutting in time and frequency.
% FORMAT D = spm_eeg_crop(S)
%
% S        - optional input struct
% (optional) fields of S:
% D        - MEEG object or filename of M/EEG mat-file with epoched data
% timewin  - time window to retain (in PST ms)
% freqwin  - frequency window to retain
% channels - cell array of channel labels or 'all'.
%
% Output:
% D        - MEEG object (also written on disk)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_crop.m 4296 2011-04-07 12:51:48Z vladimir $

SVNrev = '$Rev: 4296 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Crop M/EEG data'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

isTF = strncmpi(D.transformtype,'TF',2);

if ~isfield(S, 'timewin')
    S.timewin = spm_input('Time window (ms)', '+1', 'r', num2str(1000*[D.time(1) D.time(end)]), 2);
end

if  isTF && ~isfield(S, 'freqwin')
    S.freqwin = spm_input('Frequency window (Hz)', '+1', 'r', num2str([D.frequencies(1) D.frequencies(end)]), 2);
end

timeind = D.indsample(1e-3*(min(S.timewin))):D.indsample(1e-3*(max(S.timewin)));
if isempty(timeind) || any(isnan(timeind))
    error('Selected time window is invalid.');
end

if isTF
    freqind = D.indfrequency(min(S.freqwin)):D.indfrequency(max(S.freqwin));
    if isempty(freqind) || any(isnan(freqind))
        error('Selected frequency window is invalid.');
    end
end

if D.nchannels > 1
    if ~isfield(S, 'channels')
        [selection, ok]= listdlg('ListString', D.chanlabels, 'SelectionMode', 'multiple' ,'Name', 'Select channels' , 'ListSize', [400 300]);
        if ~ok
            return;
        end
        
        S.channels = D.chanlabels(selection);
    end
    
    if isequal(S.channels, 'all')
        chanind = 1:D.nchannels;
    else
        chanind = spm_match_str(D.chanlabels, S.channels);
    end
else
    chanind = 1;
end


%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
if isTF
    Dnew = clone(D, ['p' fnamedat(D)], [length(chanind) length(freqind) length(timeind) D.ntrials]);
    Dnew = frequencies(Dnew, [], D.frequencies(freqind));
else
    Dnew = clone(D, ['p' fnamedat(D)], [length(chanind) length(timeind) D.ntrials]);
end

Dnew = timeonset(Dnew, D.time(timeind(1)));

Dnew = chanlabels(Dnew, [], D.chanlabels(chanind));
Dnew = badchannels(Dnew, [], badchannels(D, chanind));
Dnew = chantype(Dnew, [], chantype(D, chanind));
Dnew = units(Dnew, [], units(D, chanind));
Dnew = coor2D(Dnew, [], coor2D(D, chanind));

if isequal(Dnew.type, 'continuous')
   ev = Dnew.events;
   ev = ev([ev.time]>=Dnew.time(1) & [ev.time]<=Dnew.time(end));
   Dnew   = events(Dnew, 1, ev);
end

%-Copy data
%--------------------------------------------------------------------------
spm_progress_bar('Init', D.ntrials, 'Trials copied');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end

for i = 1:D.ntrials
    
    if isTF
        Dnew(:, :, :, i) =  D(chanind, freqind, timeind, i);
    else
        Dnew(:, :, i) =  D(chanind, timeind, i);
    end
    
    if D.trialonset(i) ~= 0
        Dnew = trialonset(Dnew, i,  D.trialonset(i)+ D.time(timeind(1))-D.time(1));
    end
    
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end  %

spm_progress_bar('Clear');

%-Save the new M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','Crop M/EEG data: done'); spm('Pointer','Arrow');

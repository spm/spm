function [D, alljumps] = spm_eeg_remove_jumps(S)
% Remove "jumps" (discontinuities) from the M/EEG raw signal
% FORMAT [D, alljumps] = spm_eeg_remove_jumps(S)
%
% INPUT:
% S          - struct (optional)
% (optional) fields of S:
% D          - filename
% channels   - vector of channel numbers to be filtered, or cell array of labels
%              default = all MEEG channels.
% threshold  - threshold, default = 3000 fT (3pT)
% stdthreshold - if present overrides the threshold field and specifies the
%               threshold in terms of standard deviation
% remove     - if set to zero, jumps are detected and marked, but not
%              removed
%
% OUTPUT:
% D          - MEEG object
% alljumps   - onset of "jumps" in seconds, summarized across channels into
%              timebins of 0.1 s length
%__________________________________________________________________________
%
% This function removes "jumps" (discontinuities) from the EEG/MEG raw
% signal, based on an absolute threshold, and filters the signal derivative
% over 20 timepoints.
% Such jumps occur with squid resetting and when acquisition is stopped
% with the "abort" button.
% This procedure is necessary before performing highpass filtering on the
% continuous data.
% Timestamps for the jumps are returned and at the same time recorded as
% event markers. Epochs containing such jumps should be rejected as they are
% affected by ringing from analogue filters in the recording system.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Dominik R Bach
% $Id: spm_eeg_remove_jumps.m 3463 2009-10-14 11:30:05Z vladimir $

% Input parameters
%==========================================================================

% get file
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select M/EEG mat file');
    S.D = D;
end

if ischar(D)
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end


% get threshold
%--------------------------------------------------------------------------
try
    thr      = S.threshold;
catch
    thr      = 3e-12;       % 3000 fT
end

% get channels
%--------------------------------------------------------------------------
try
    channels = S.channels;
    if ischar(channels)
        channels = {channels};
    end
    if iscell(channels)
        channels = spm_match_str(D.chanlabels, channel);
    end
catch
    channels = meegchannels(D);
end

% remove or detect? 
% -------------------------------------------------------------------------
try
    remove = S.remove;
    if ~isnumeric(remove) || ~ismember(remove, 1:2)
        remove = 1;
    end;
catch
    remove = 1;
end;


% Detect and remove jumps
%==========================================================================

% progress report
%--------------------------------------------------------------------------
[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'M/EEG remove jumps', 0);
spm('Pointer', 'Watch');

% generate new meeg object with new filenames
%--------------------------------------------------------------------------
S1   = [];
S1.D = D;
S1.newname = ['j' fname(D)];

Dnew = spm_eeg_copy(S1);


if isequal(D.type, 'continuous')
    % determine block size, depending on available memory
    %--------------------------------------------------------------------------
    try
        % 2/3 of largest block of contiguous memory, for Windows platforms
        evalc('memsz=2/3*feature(''memstats'');');
    catch
        % 20 MB for all other platforms
        memsz = 20*1024*1024;
    end
    datasz    = numel(channels)*nsamples(D)*8;           % datapoints x 8 bytes
    blknum    = ceil(datasz/memsz);
    blksz     = ceil(numel(channels)/blknum);
    spm_progress_bar('Init', numel(channels), 'Channels filtered');
else
    blknum    = D.ntrials;
    spm_progress_bar('Init', blknum, 'Trials filtered');
end

% filter blocks of channels
%--------------------------------------------------------------------------
chncnt = 1;
for blk = 1:blknum

    if isequal(D.type, 'continuous')
        % load original data blockwise into workspace
        %----------------------------------------------------------------------
        blkchan = chncnt:(min(numel(channels), chncnt+blksz-1));
        Dtemp   = D(channels(blkchan),:,1);
        chncnt  = chncnt + blksz;
    else
        blkchan = 1:length(channels);
        Dtemp   = D(channels,:,blk);
    end

    jumps_fixed = 0;
    % loop through channels within blocks
    %----------------------------------------------------------------------
    for ch = 1:numel(blkchan)

        % find jumps in derivative
        %------------------------------------------------------------------
        dat   = Dtemp(ch,:,1);
        ddat  = diff(dat);
        if ~isfield(S, 'stdthreshold')
            jumps = find(abs(ddat) > thr);
        else
            jumps = find(abs(ddat) > S.stdthreshold*std(dat));
        end
        if ~isempty(jumps) && remove

            % collapse jumps than are closer than 10 timepoints apart
            %--------------------------------------------------------------
            if numel(jumps)>2, jumps(find(diff(jumps) < 10) + 1) = []; end;

            % set derivative to median over +- 20 timepoints
            %--------------------------------------------------------------
            for i = 1:numel(jumps)
                ind = jumps(i) + (-20:20);
                if ind(1) < 1
                    ind = ind + 1 - ind(1);
                end
                if ind(end) > length(ddat)
                    ind = ind - (ind(end) - length(ddat));
                end
                ddat(ind) = median(ddat(ind));
            end

            % reconstruct data
            %--------------------------------------------------------------
            dat(2:end) = cumsum(ddat) + dat(1);

            jumps_fixed       = 1;
        end
        
        % store jumps onsets and filtered data
        %------------------------------------------------------------------
        if isequal(D.type, 'continuous')
            jmps{blkchan(ch), 1} = jumps;
        else
            jmps{blkchan(ch), blk} = jumps;
        end

        Dtemp(ch,:)       = dat;

        % indicate progress
        %--------------------------------------------------------------
        if isequal(D.type, 'continuous')
            spm_progress_bar('Set', blkchan(ch));
        end
    end

    if jumps_fixed
        % write filtered data blockwise in new data file
        %----------------------------------------------------------------------
        if isequal(D.type, 'continuous')
            Dnew(channels(blkchan),:,1) = Dtemp;
        else
            Dnew(channels,:,blk) = Dtemp;

        end
    end

    if ~isequal(D.type, 'continuous')
        spm_progress_bar('Set', blk);
    end
end

spm_progress_bar('Clear');


% Insert artefact timepoints as event markers of type "artefact"
%==========================================================================
alljumps = cell(1, D.ntrials);

for n = 1:D.ntrials
    % summarize jumps across channels into 0.1 s timebins, expressed in seconds
    %--------------------------------------------------------------------------
    alljumps{n} = unique( ceil(sort(cell2mat(jmps(:, n)')) / fsample(D) * 10) ) / 10;

    
    trialonset  = D.trialonset(n);
    if iscell(trialonset)
        trialonset = trialonset{1};
    end
    
    if isempty(trialonset)
        trialonset = D.timeonset;
    end
    

    ev = events(Dnew, n);

    if iscell(ev)
        ev = ev{1};
    end

    Nevents = numel(ev);
    for i=1:numel(alljumps{n})
        ev(Nevents+i).type     = 'artefact';
        ev(Nevents+i).value    = 'jump';
        ev(Nevents+i).duration = [];
        ev(Nevents+i).time     = alljumps{n}(i)+trialonset;
    end
    
    if ~isempty(ev)
        [tevent, I] = sort([ev.time]);
        ev = ev(I);
        Dnew = events(Dnew, n, ev);
    end
end

if numel(alljumps) == 1
    alljumps = alljumps{1};
end

% Save new meeg object
%==========================================================================
D = Dnew;
D = D.history('spm_eeg_remove_jumps', S);
save(D);

% report
%--------------------------------------------------------------------------
fprintf('%s: Found %d jumps.\n', fname(D), length(spm_vec(alljumps)));            %-#

spm('Pointer', 'Arrow');
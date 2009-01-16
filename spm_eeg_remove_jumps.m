function [D, alljumps] = spm_eeg_remove_jumps(S)
% Remove "jumps" (discontinuities) from the M/EEG raw signal 
% FORMAT [D, alljumps] = spm_eeg_remove_jumps(S)
%
% INPUT:
% S          - struct (optional)
% (optional) fields of S:
% D          - filename
% channels   - vector of channel numbers to be filtered, default = all
% threshold  - threshold, default = 3000 fT (3pT)
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
% $Id: spm_eeg_remove_jumps.m 2610 2009-01-16 14:54:13Z guillaume $

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

if ~strcmpi(type(D), 'continuous')
    error('This function is meant to be performed on continuous data only.');
end;

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
catch
    channels = 1:nchannels(D);
end

% Detect and remove jumps
%==========================================================================

% progress report
%--------------------------------------------------------------------------
[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'M/EEG remove jumps', 0);
spm('Pointer', 'Watch');
spm_progress_bar('Init', numel(channels), 'Channels filtered');

% generate new meeg object with new filenames
%--------------------------------------------------------------------------
Dnew = clone(D, ['j' fnamedat(D)], [D.nchannels D.nsamples D.ntrials]);

% copy channels not to be filtered
%--------------------------------------------------------------------------
ind = setdiff(1:D.nchannels, channels);
if ~isempty(ind)
    for i = 1:length(ind)
        Dnew(ind(i),:,1) = D(ind(i),:,1);
    end
end

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

% filter blocks of channels
%--------------------------------------------------------------------------
chncnt = 1;
for blk = 1:blknum
    
    % load original data blockwise into workspace 
    %----------------------------------------------------------------------
    blkchan = chncnt:(min(numel(channels), chncnt+blksz-1));
    Dtemp   = D(channels(blkchan),:,1);
    chncnt  = chncnt + blksz;
    
    % loop through channels within blocks
    %----------------------------------------------------------------------
    for ch = 1:numel(blkchan)
        
        % find jumps in derivative
        %------------------------------------------------------------------
        dat   = Dtemp(ch,:,1);
        ddat  = diff(dat);
        jumps = find(abs(ddat) > thr);
        if ~isempty(jumps)
            
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
        end
        
        % store jumps onsets and filtered data
        %------------------------------------------------------------------
        jmps{blkchan(ch)} = jumps;
        Dtemp(ch,:)       = dat;
        
        % indicate progress
        %------------------------------------------------------------------
        spm_progress_bar('Set', blkchan(ch)); 
    end

    % write filtered data blockwise in new data file
    %----------------------------------------------------------------------
    Dnew(channels(blkchan),:,1) = Dtemp;
end

spm_progress_bar('Clear');

% summarize jumps across channels into 0.1 s timebins, expressed in seconds
%--------------------------------------------------------------------------
alljumps = unique( ceil(sort(cell2mat(jmps)) / fsample(D) * 10) ) / 10;

% Insert artefact timepoints as event markers of type "artefact"
%==========================================================================
ev = events(Dnew, 1);

Nevents = numel(ev);
for i=1:numel(alljumps)
    ev(Nevents+i).type     = 'artefact';
    ev(Nevents+i).value    = [];
    ev(Nevents+i).duration = [];
    ev(Nevents+i).time     = alljumps(i);
end
[tevent, I] = sort([ev.time]);
ev = ev(I);

Dnew = events(Dnew, 1, ev);

% Save new meeg object
%==========================================================================
D = Dnew;
D = D.history('spm_eeg_remove_jumps', S);
save(D);

% report
%--------------------------------------------------------------------------
fprintf('%s: Found %d jumps.\n', fname(D), numel(alljumps));            %-#

spm('Pointer', 'Arrow');
function D = spm_eeg_remove_jumps(S)
% Remove "jumps" (discontinuities) from the M/EEG raw signal
% FORMAT [D, alljumps] = spm_eeg_remove_jumps(S)
%
% INPUT:
% S          - struct (optional)
% fields of S:
% D          - filename
% channels          - cell array of channel names. Can include generic
%                         wildcards: 'All', 'EEG', 'MEG' etc.
% threshold  - threshold, default = 3000 fT (3pT)
% stdthreshold - if present overrides the threshold field and specifies the
%               threshold in terms of standard deviation
% prefix       - prefix for the output dataset (default - 'j')
%
% OUTPUT:
% D          - MEEG object
%__________________________________________________________________________
%
% This function removes "jumps" (discontinuities) from the EEG/MEG raw
% signal, based on an absolute threshold, and filters the signal derivative
% over 20 timepoints.
% Such jumps occur with squid resetting and when acquisition is stopped
% with the "abort" button.
% This procedure is necessary before performing highpass filtering on the
% continuous data.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Dominik R Bach
% $Id: spm_eeg_remove_jumps.m 5865 2014-01-31 11:48:30Z vladimir $

% Input parameters
%==========================================================================

SVNrev = '$Rev: 5865 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Remove jumps'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix      = 'j';           end
if ~isfield(S, 'threshold'),    S.threshold   = 3000;          end

D = spm_eeg_load(S.D);

channels = D.selectchannels(S.channels);

% Detect and remove jumps
%==========================================================================
spm('Pointer', 'Watch');


Dnew = copy(D, [S.prefix fname(D)]);


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

win = round(0.05*D.fsample);

% filter blocks of channels
%--------------------------------------------------------------------------
chncnt = 1;
alljumps = 0;
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
            jumps = find(abs(ddat) > S.threshold);
        else
            jumps = find(abs(ddat) > S.stdthreshold*std(ddat));
        end
      
        while  ~isempty(jumps)
%%
            % collapse jumps than are closer than 10 timepoints apart
            %--------------------------------------------------------------
            if numel(jumps)>2, jumps(find(diff(jumps) < 10) + 1) = []; end;

            % set derivative to median over +- 20 timepoints
            %--------------------------------------------------------------
            for i = 1:numel(jumps)
                ind = jumps(i) + (-win:win);
                if ind(1) < 1
                    ind = ind + 1 - ind(1);
                end
                if ind(end) > length(ddat)
                    ind = ind - (ind(end) - length(ddat));
                end
                
                if abs(median(ddat(ind)))<max(abs(ddat(ind)))
                    ddat(ind) = median(ddat(ind));
                else %Just make it flat
                    ddat(ind) = 0;
                end
                
                jumps_fixed       = jumps_fixed + 1;
            end

            % reconstruct data
            %--------------------------------------------------------------
            dat(2:end) = cumsum(ddat) + dat(1);      
            
            if ~isfield(S, 'stdthreshold')
                jumps = find(abs(ddat) > S.threshold);
            else
                jumps = find(abs(ddat) > S.stdthreshold*std(ddat));
            end
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

    alljumps = alljumps + jumps_fixed;
    
    if ~isequal(D.type, 'continuous')
        spm_progress_bar('Set', blk);
    end
end

spm_progress_bar('Clear');


% Save new meeg object
%==========================================================================
D = Dnew;
D = D.history('spm_eeg_remove_jumps', S);
save(D);

% report
%--------------------------------------------------------------------------
fprintf('%s: Found %d jumps.\n', fname(D), alljumps);            %-#

spm('Pointer', 'Arrow');
function D = spm_eeg_cont_power(S)
% Compute power of continuous M/EEG data
% FORMAT D = spm_eeg_cont_power(S)
%
% This function computes power from band-pass filtered data using hilbert
% transform. Can also be used as a template fof any kind of computation on
% continuous data.
%
% S           - input structure (optional)
% (optional) fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file
%
% D           - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_cont_power.m 4286 2011-04-04 13:35:04Z vladimir $

SVNrev = '$Rev: 4286 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG power'); spm('Pointer', 'Watch');

if nargin == 0
    S = [];
end

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

%
%--------------------------------------------------------------------------

% generate new meeg object with new filenames
Dnew = clone(D, ['h' fnamedat(D)], [D.nchannels D.nsamples D.ntrials]);

% determine channels to process
chanind = D.meegchannels;


% continuous data
spm_progress_bar('Init', nchannels(D), 'Channels processed'); drawnow;
if nchannels(D) > 100, Ibar = floor(linspace(1, nchannels(D),100));
else Ibar = [1:nchannels(D)]; end

% work on blocks of channels
% determine blocksize
% determine block size, dependent on memory
memsz  = spm('Memory');
datasz = nchannels(D)*nsamples(D)*8; % datapoints x 8 bytes per double value
blknum = ceil(datasz/memsz);
blksz  = ceil(nchannels(D)/blknum);
blknum = ceil(nchannels(D)/blksz);

% now filter blocks of channels
chncnt=1;
for blk=1:blknum
    % load old meeg object blockwise into workspace
    blkchan=chncnt:(min(nchannels(D), chncnt+blksz-1));
    if isempty(blkchan), break, end
    tempdata=D(blkchan,:,1);
    chncnt=chncnt+blksz;
    %loop through channels
    for j = 1:numel(blkchan)
        
        if ismember(blkchan(j), chanind)
            % ********** this is actually the processing part
            tempdata(j, :) = abs(spm_hilbert(tempdata(j,:)));
            % **********
        end
        
        if ismember(j, Ibar), spm_progress_bar('Set', blkchan(j)); end
        
    end
    
    % write tempdata to Dnew
    Dnew(blkchan,:,1)=tempdata;
    clear tempdata;   
end;


spm_progress_bar('Clear');

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = Dnew;
D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG power: done'); spm('Pointer', 'Arrow');


function D = spm_eeg_reduce(S)
% Apply data reduction to M/EEG dataset
% FORMAT D = spm_eeg_reduce(S)
% S                     - input structure
%
% fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file with
%   
%   S.channels          - cell array of channel names. Can include generic
%                         wildcards: 'All', 'EEG', 'MEG' etc
%   S.method           - name for the spectral estimation to use. This
%                        corresponds to the name of a plug-in function that comes
%                        after 'spm_eeg_reduce_' prefix.
%   S.settings         - plug-in specific settings
%
% Output:
% D                     - M/EEG object 
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_reduce.m 4911 2012-09-07 15:21:40Z vladimir $

SVNrev = '$Rev: 4911 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG reduce'); spm('Pointer','Watch');

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


%-Configure the analysis
%--------------------------------------------------------------------------
if ~isfield(S, 'channels')
    S.channels = 'all';
end

badind = D.badchannels;

chanind = setdiff(D.selectchannels(S.channels), badind);

if isempty(chanind)
    error('No channels selected.');
end

if ~isfield(S, 'method')
    S.method = 'pca';
    S.settings.ncomp = min(length(floor(chanind)/2), 100);
end

if isfield(S, 'settings')
    S1 = S.settings;
else
    S1 = [];
end

S1.D       = D;
S1.chanind = chanind;

montage = feval(['spm_eeg_reduce_' S.method], S1);

% This is to discard bad channels but keep other channels (like non MEEG).

if ~isempty(badind)
    montage.labelorg = [montage.labelorg(:) D.chanlabels(badind)'];
    montage.tra(end, end+length(badind)) = 0;
end


S1 = [];
S1.D = D;
S1.montage = montage;
S1.keepothers = 'yes';
S1.updatehistory  = 0;
D = spm_eeg_montage(S1);


%-Save new M/EEG dataset(s)
%--------------------------------------------------------------------------
D = D.history('spm_eeg_reduce', S);
save(D);
%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG reduce: done'); spm('Pointer','Arrow');

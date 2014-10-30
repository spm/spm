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
%   S.conditions          - cell array of condition trial names. 
%   S.method           - name for the spectral estimation to use. This
%                        corresponds to the name of a plug-in function that comes
%                        after 'spm_eeg_reduce_' prefix.
%   S.keepothers       - keep the other (not involved) channels
%   S.settings         - plug-in specific settings
%   S.timewin          - time windows or interest
%   S.prefix           - prefix for the output file (default - 'R')
%
% Output:
% D                     - M/EEG object 
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_reduce.m 6252 2014-10-30 16:06:30Z vladimir $

SVNrev = '$Rev: 6252 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG reduce'); spm('Pointer','Watch');


%-Configure the analysis
%--------------------------------------------------------------------------
if ~isfield(S, 'channels'),   S.channels = 'all';             end
if ~isfield(S, 'conditions'), S.conditions.all = 1;           end
if ~isfield(S, 'timewin'),    S.timewin  = [-Inf Inf];        end
if ~isfield(S, 'keepothers'), S.keepothers  = true;           end
if ~isfield(S, 'prefix'),     S.prefix   = 'R';               end


D = spm_eeg_load(S.D);

badind = D.badchannels;

chanind = setdiff(D.selectchannels(S.channels), badind);

% if isempty(chanind)
%     error('No channels selected.');
% end

%%%%%%%%%
% MWW
samples = {};
for i = 1:size(S.timewin, 1)
    samples{i} = D.indsample(S.timewin(i, 1)):D.indsample(S.timewin(i, 2));
end

if isfield(S.conditions, 'all')
    trials = 1:D.ntrials;
else    
    trials = D.indtrial(S.conditions, 'GOOD');
    if isempty(trials)
        error('No trials matched the selection, check the specified condition labels');
    end
end
%%%%%%%%%%


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
S1.trials = trials; 
S1.samples = samples; 
montage = feval(['spm_eeg_reduce_' S.method], S1);

% This is to discard bad channels but keep other channels (like non MEEG).

if ~isempty(badind)
    montage.labelorg = [montage.labelorg(:); D.chanlabels(badind)']; % added semicolon - MWW
    montage.tra(end, end+length(badind)) = 0;
    
    if isfield(montage, 'chantypeorg')
        montage.chantypeorg = [montage.chantypeorg(:); lower(D.chantype(badind))'];
    end
    
    if isfield(montage, 'chanunitorg')
        montage.chanunitorg = [montage.chanunitorg(:); D.units(badind)'];
    end
end

S1 = [];
S1.D = D;
S1.montage = montage;
S1.keepothers = S.keepothers; 
S1.prefix = S.prefix;
S1.updatehistory  = 0;
D = spm_eeg_montage(S1);

D = chantype(D, 1:length(montage.labelnew), montage.chantypenew);

%-Save new M/EEG dataset(s)
%--------------------------------------------------------------------------
D = D.history('spm_eeg_reduce', S);
save(D);
%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG reduce: done'); spm('Pointer','Arrow');

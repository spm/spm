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
%   S.keeporig         - keep the original unreduced channels (1) or remove
%                        (0, default).
%   S.keepothers       - keep the other (not involved) channels
%   S.settings         - plug-in specific settings
%   S.timewin          - time windows or interest
%   S.prefix           - prefix for the output file (default - 'R')
%
% Output:
% D                     - M/EEG object 
%__________________________________________________________________________
% Copyright (C) 2012-2016 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_reduce.m 6843 2016-07-28 10:55:47Z vladimir $

SVNrev = '$Rev: 6843 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG reduce'); spm('Pointer','Watch');


%-Configure the analysis
%--------------------------------------------------------------------------
if ~isfield(S, 'channels'),   S.channels = 'all';             end
if ~isfield(S, 'conditions'), S.conditions.all = 1;           end
if ~isfield(S, 'timewin'),    S.timewin  = [-Inf Inf];        end
if ~isfield(S, 'keeporig'),   S.keeporig = false;          end
if ~isfield(S, 'keepothers'), S.keepothers  = true;           end
if ~isfield(S, 'prefix'),     S.prefix   = 'R';               end


if iscell(S.D)
    badind = [];
    for i = 1:numel(S.D)
        DD{i} = spm_eeg_load(S.D{i});
        badind = [badind DD{i}.badchannels];    
    end
    D = DD{1};
    badind = unique(badind);
else
    D = spm_eeg_load(S.D);
    DD{1}  = D;
    badind = D.badchannels;
end

chanind = setdiff(D.selectchannels(S.channels), badind);

if isempty(chanind)
    error('No channels selected.');
end

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

if numel(DD)>1
    S1.D = DD;
else
    S1.D = D;
end

S1.chanind = chanind; 
S1.trials = trials; 
S1.samples = samples; 
montage = feval(['spm_eeg_reduce_' S.method], S1);

if ~isempty(badind)
    montage.labelorg = [montage.labelorg(:); D.chanlabels(badind)'];
    if S.keeporig
        montage.tra((end+1):(end+length(badind)), (end+1):(end+length(badind))) = eye(length(badind));     
    else
        montage.tra(end, end+length(badind)) = 0;
    end
    
    if isfield(montage, 'chantypeorg')
        montage.chantypeorg = [montage.chantypeorg(:); lower(D.chantype(badind))'];
    end
    
    if isfield(montage, 'chanunitorg')
        montage.chanunitorg = [montage.chanunitorg(:); D.units(badind)'];
    end
end

% This is to discard bad channels but keep other channels (like non MEEG).
if S.keeporig
    montage.labelnew = [montage.labelnew(:); D.chanlabels(chanind)'];
    montage.tra((end+1):(end+length(chanind)), 1:length(chanind)) = eye(length(chanind));    
    
    if isfield(montage, 'chantypenew')
        montage.chantypenew = [montage.chantypenew(:); lower(D.chantype([badind chanind]))'];
    end
    
    if isfield(montage, 'chanunitnew')
        montage.chanunitnew = [montage.chanunitnew(:); D.units([badind chanind])'];
    end
end

% Reorder as in the original file. This might be handy when reducing and
% then grand-averaging or merging across subjects
[sel1, sel2] = spm_match_str(D.chanlabels, montage.labelnew);
sortind = [spm_vec(sel2); spm_vec(setdiff(1:length(montage.labelnew), sel2))];

montage.labelnew = montage.labelnew(sortind);
montage.tra = montage.tra(sortind, :);
if isfield(montage, 'chantypenew')
    montage.chantypenew = montage.chantypenew(sortind);
end
if isfield(montage, 'chanunitnew')
    montage.chanunitnew = montage.chanunitnew(sortind);
end

S1 = [];
S1.montage = montage;
S1.keepothers = S.keepothers; 
S1.prefix = S.prefix;
S1.updatehistory  = 0;
for i = 1:numel(DD)
    S1.D = DD{i};
    D = spm_eeg_montage(S1);
    
    %-Save new M/EEG dataset(s)
    %--------------------------------------------------------------------------
    D = D.history('spm_eeg_reduce', S);
    save(D);
    
    DD{i} = D;
end

if numel(DD)>1
    D = DD;
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG reduce: done'); spm('Pointer','Arrow');

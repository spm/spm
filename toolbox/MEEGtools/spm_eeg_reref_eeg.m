function [D, S] = spm_eeg_reref_eeg(S)
% Rereference EEG data to new reference channel(s)
% FORMAT [D, S] = spm_eeg_reref_eeg(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.refchan          - New reference channel indices or labels
%                        ('average' can be used as shortcut)
%
% Output:
% D                    - MEEG object (also written on disk)
% S                    - record of parameters, including montage
%__________________________________________________________________________
%
% spm_eeg_reref_eeg re-references any EEG data within an MEEG dataset, by
% calling spm_eeg_montage with appropriate montage, excluding bad channels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Rik Henson (updated so that bad channels included in tra but not used)
% $Id: spm_eeg_reref_eeg.m 4119 2010-11-11 12:12:00Z vladimir $

SVNrev = '$Rev: 4119 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','EEG reref'); spm('Pointer','Watch');

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

% Get indices for just EEG channels and remove any bad channels
%--------------------------------------------------------------------------
eegchan  = D.meegchannels('EEG');
goodchan = setdiff(eegchan, D.badchannels);

if isempty(goodchan)
    error('No EEG channels that are not marked as bad...!?')
end

%-Get reference channel indices 
%--------------------------------------------------------------------------
if ~isfield(S, 'refchan')
      [selection, ok]= listdlg('ListString', D.chanlabels(goodchan), 'SelectionMode', 'multiple' ,'Name', 'Select reference channels' , 'ListSize', [400 300]);
    if ~ok
        return;
    end
    
    if length(selection) == length(goodchan)
        S.refchan = 'average';
    else
        S.refchan = D.chanlabels(goodchan(selection));
    end
end

if strcmpi(S.refchan,'average')
    refchan = goodchan;
elseif iscell(S.refchan)
    refchan = D.indchannel(S.refchan);
elseif isnumeric(S.refchan)
    refchan = S.refchan;
end

refind  = find(ismember(eegchan, refchan));
goodind = find(ismember(eegchan, goodchan));
badind  = find(ismember(eegchan, D.badchannels));

if length(refind) ~= length(refchan)
    error('Not all S.refchan are valid indices for (non-bad) EEG channels')
end

tra                 = eye(length(eegchan));
tra(goodind,refind) = tra(goodind,refind) - 1/length(refchan);
tra(badind,refind)  = tra(badind,refind)  - 1/length(refchan);

S1=[];
S1.D = D;
S1.montage.labelorg = D.chanlabels(eegchan);
S1.montage.labelnew = D.chanlabels(eegchan);
S1.montage.tra = tra;
S1.keepothers = 'yes';
S1.updatehistory = 0;
D = spm_eeg_montage(S1);

%-Update history (not necessary; leave as call to spm_eeg_montage?) Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = D.history('spm_eeg_reref_eeg', S);

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','EEG re-reference: done'); spm('Pointer','Arrow');

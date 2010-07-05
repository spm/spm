function [D, montage] = spm_eeg_montage(S)
% Rereference EEG data to new reference channel(s)
% FORMAT [D, montage] = spm_eeg_montage(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.montage          - (optional)
%     A montage is specified as a structure with the fields:
%     montage.tra      - MxN matrix
%     montage.labelnew - Mx1 cell-array: new labels
%     montage.labelorg - Nx1 cell-array: original labels
%     or as a filename of a mat-file containing the montage structure.
%   S.keepothers       - 'yes'/'no': keep or discart the channels not
%                        involved in the montage [default: 'yes']
%   S.blocksize        - size of blocks used internally to split large
%                        continuous files [default ~20Mb]
%   S.updatehistory    - if 0 the history is not updated (useful for
%                        functions that use montage functionality.
%
% Output:
% D                    - MEEG object (also written on disk)
% montage              - the applied montage
%__________________________________________________________________________
%
% spm_eeg_montage applies montage provided or specified by the user to
% data and sensors of an MEEG file and produces a new file. It works
% correctly when no sensors are specified or when data and sensors are
% consistent (which is ensured by spm_eeg_prep_ui).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Robert Oostenveld, Stefan Kiebel
% $Id: spm_eeg_montage.m 3970 2010-07-05 17:26:16Z vladimir $

SVNrev = '$Rev: 3970 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG montage'); spm('Pointer','Watch');

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

if strncmp(D.transformtype, 'TF', 2)
    error('Montage cannot be applied to time-frequency data');
end

%-Get size of blocks used internally to split large continuous files
%--------------------------------------------------------------------------
try
    S.blocksize;
catch
    S.blocksize = 655360; %20 Mb
end

%-Get montage
%--------------------------------------------------------------------------
if ~isfield(S, 'montage')
    res = spm_input('How to specify the montage?', '+1', 'gui|file');
    switch res
        case 'gui'
            montage.labelorg = D.chanlabels;
            montage.labelnew = montage.labelorg;
            montage.tra      = eye(length(montage.labelorg));
            montage.labelnew = montage.labelorg;
            S.montage        = spm_eeg_montage_ui(montage);
        case 'file'
            [S.montage, sts] = spm_select(1, 'mat', 'Select a montage file');
            if ~sts, montage = []; return; end
    end
end

if ischar(S.montage)
    try
        montage = load(S.montage);
    catch
        error(sprintf('Could not read montage file %s.',S.montage));
    end
    if ismember('montage', fieldnames(montage))
        S.montage = montage.montage;
    else
        error(sprintf('Invalid montage file %s.',S.montage));
    end
end

montage = S.montage;

if ~all(isfield(montage, {'labelnew', 'labelorg', 'tra'})) || ...
        any([isempty(montage.labelnew), isempty(montage.labelorg), isempty(montage.tra)]) || ...
        length(montage.labelnew) ~= size(montage.tra, 1) || length(montage.labelorg) ~= size(montage.tra, 2)
    error('Insufficient or illegal montage specification.');
end

% select and discard the columns that are empty
selempty         = find(~all(montage.tra == 0, 1));
montage.tra      = montage.tra(:, selempty);
montage.labelorg = montage.labelorg(selempty);

% add columns for the channels that are not involved in the montage
[add, ind]       = setdiff(D.chanlabels, montage.labelorg);
chlab            = D.chanlabels;
add              = chlab(sort(ind));

%-Get keepothers input field
%--------------------------------------------------------------------------
if ~isfield(S, 'keepothers')
    if ~isempty(add)
        S.keepothers  = spm_input('Keep the other channels?', '+1', 'yes|no');
    else
        S.keepothers = 'yes';
    end
end

m = size(montage.tra, 1);
n = size(montage.tra, 2);
k = length(add);
if strcmp(S.keepothers, 'yes')
    montage.tra((m+(1:k)), (n+(1:k))) = eye(k);
    montage.labelnew = cat(1, montage.labelnew(:), add(:));
else
    montage.tra = [montage.tra zeros(m, k)];
end
montage.labelorg = cat(1, montage.labelorg(:), add(:));

% determine whether all channels are unique
m = size(montage.tra,1);
n = size(montage.tra,2);
if length(unique(montage.labelnew))~=m
    error('not all output channels of the montage are unique');
end
if length(unique(montage.labelorg))~=n
    error('not all input channels of the montage are unique');
end

% determine whether all channels that have to be rereferenced are available
if length(intersect(D.chanlabels, montage.labelorg)) ~= n
    error('not all channels that are used in the montage are available');
end

% reorder the columns of the montage matrix
[selchan, selmont] = spm_match_str(D.chanlabels, montage.labelorg);
montage.tra        = montage.tra(:,selmont);
montage.labelorg   = montage.labelorg(selmont);

spm('Pointer', 'Watch');

%-Generate new MEEG object with new filenames
%--------------------------------------------------------------------------
Dnew = clone(D, ['M' fnamedat(D)], [m D.nsamples D.ntrials], 1);

nblocksamples = floor(S.blocksize/max(D.nchannels, m));

if D.nsamples <= nblocksamples
    nblocks = 1;
    nblocksamples = D.nsamples;
else
    nblocks = floor(D.nsamples./nblocksamples);
end

if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
elseif  D.ntrials == 1, Ibar = [1:ceil(D.nsamples./nblocksamples)];
else Ibar = [1:D.ntrials]; end

spm_progress_bar('Init', Ibar(end), 'applying montage');

for i = 1:D.ntrials
    for j = 1:nblocks
        
        Dnew(:, ((j-1)*nblocksamples +1) : (j*nblocksamples), i) = ...
            montage.tra * squeeze(D(:, ((j-1)*nblocksamples +1) : (j*nblocksamples), i));
        
        if D.ntrials == 1 && ismember(j, Ibar)
            spm_progress_bar('Set', j);
        end
    end
    
    if mod(D.nsamples, nblocksamples) ~= 0
        Dnew(:, (nblocks*nblocksamples +1) : D.nsamples, i) = ...
            montage.tra * squeeze(D(:, (nblocks*nblocksamples +1) : D.nsamples, i));
    end
    
    
    if D.ntrials>1 && ismember(i, Ibar)
        spm_progress_bar('Set', i);
    end
end

Dnew = chanlabels(Dnew, 1:Dnew.nchannels, montage.labelnew);

% Transfer bad flags in case there are channels with the
% same name in both files.
if  ~isempty(D.badchannels)
    sel = spm_match_str(Dnew.chanlabels, D.chanlabels(D.badchannels));
    if ~isempty(sel)
        Dnew = badchannels(Dnew, sel, 1);
    end
end

% Set channel types to default
S1 = [];
S1.task = 'defaulttype';
S1.D = Dnew;
S1.updatehistory = 0;
Dnew = spm_eeg_prep(S1);

%-Apply montage to sensors
%--------------------------------------------------------------------------
sensortypes = {'MEG', 'EEG'};
for i = 1:length(sensortypes)
    sens = D.sensors(sensortypes{i});
    if ~isempty(sens) && length(intersect(sens.label, montage.labelorg))==length(sens.label) &&...
            ~isequal(sensortypes{i}, 'EEG') % EEG disabled for now
        sensmontage = montage;
        sel1 = spm_match_str(sensmontage.labelorg, sens.label);
        sensmontage.labelorg = sensmontage.labelorg(sel1);
        sensmontage.tra = sensmontage.tra(:, sel1);
        selempty  = find(all(sensmontage.tra == 0, 2));
        sensmontage.tra(selempty, :) = [];
        sensmontage.labelnew(selempty) = [];
        sens = ft_apply_montage(sens, sensmontage, 'keepunused', S.keepothers);
        
        if isfield(sens, 'balance') && ~isequal(sens.balance.current, 'none')
            balance = ft_apply_montage(getfield(sens.balance, sens.balance.current), sensmontage, 'keepunused', S.keepothers);
        else
            balance = sensmontage;
        end
        
        sens.balance.custom = balance;
        sens.balance.current = 'custom';
    end
    
    if ~isempty(sens) && ~isempty(intersect(sens.label, Dnew.chanlabels))
        Dnew = sensors(Dnew, sensortypes{i}, sens);
    else
        Dnew = sensors(Dnew, sensortypes{i}, []);
    end
end

% Remove any inverse solutions
if isfield(Dnew, 'inv')
    Dnew = rmfield(Dnew, 'inv');
end

%-Assign default EEG sensor positions if no positions are present or if
% default locations had been assigned before but no longer cover all the
% EEG channels.
%--------------------------------------------------------------------------
if ~isempty(Dnew.meegchannels('EEG')) && (isempty(Dnew.sensors('EEG')) ||...
        (all(ismember({'spmnas', 'spmlpa', 'spmrpa'}, Dnew.fiducials.fid.label)) && ...
        ~isempty(setdiff(Dnew.chanlabels(Dnew.meegchannels('EEG')), getfield(Dnew.sensors('EEG'), 'label')))))
    S1 = [];
    S1.task = 'defaulteegsens';
    S1.updatehistory = 0;
    S1.D = Dnew;
    
    Dnew = spm_eeg_prep(S1);
end

%-Create 2D positions for EEG by projecting the 3D positions to 2D
%--------------------------------------------------------------------------
if ~isempty(Dnew.meegchannels('EEG')) && ~isempty(Dnew.sensors('EEG'))
    S1 = [];
    S1.task = 'project3D';
    S1.modality = 'EEG';
    S1.updatehistory = 0;
    S1.D = Dnew;
    
    Dnew = spm_eeg_prep(S1);
end

%-Create 2D positions for MEG  by projecting the 3D positions to 2D
%--------------------------------------------------------------------------
if ~isempty(Dnew.meegchannels('MEG')) && ~isempty(Dnew.sensors('MEG'))
    S1 = [];
    S1.task = 'project3D';
    S1.modality = 'MEG';
    S1.updatehistory = 0;
    S1.D = Dnew;
    
    Dnew = spm_eeg_prep(S1);
end

%-Transfer the properties of channels not involved in the montage
%--------------------------------------------------------------------------
if ~isempty(add) && strcmp(S.keepothers, 'yes')
    old_add_ind = D.indchannel(add);
    new_add_ind = Dnew.indchannel(add);
    
    Dnew = badchannels(Dnew, new_add_ind, badchannels(D, old_add_ind));
    Dnew = chantype(Dnew, new_add_ind, chantype(D, old_add_ind));
end

[ok, Dnew] = check(Dnew, 'basic');

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = Dnew;

if ~isfield(S, 'updatehistory') || S.updatehistory
    D = D.history('spm_eeg_montage', S);
end

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG montage: done'); spm('Pointer','Arrow');

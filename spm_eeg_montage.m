function [D, montage] = spm_eeg_montage(S)
% Rereference EEG data to new reference channel(s)
% FORMAT [D, montage] = spm_eeg_montage(S)
%
% S                    - input structure
%  fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.mode             - 'write' (default) apply montage and create a new
%                                dataset
%                        'switch' apply online montage without writing out
%                                 the data
%                        'add'    add the montage to the set of online
%                                 montages without applying it
%                        'clear'  clear all online montages and revert to
%                                 the original state
%   S.montage          -
%     A montage is specified as a structure with the fields:
%     montage.tra      - MxN matrix
%     montage.labelnew - Mx1 cell-array: new labels
%     montage.labelorg - Nx1 cell-array: original labels
%     montage.name     - (optional) montage name when using online montage
%     or as a filename of a mat-file containing the montage structure,
%     or as the index of the online montage to use
%     or as the name of the online montage
%   S.keepothers       - keep (1) or discard (0) the channels not
%                        involved in the montage [default: 1]
%                        ignored when switching between online montages
%   S.keepsensors      - keep (1) or discard (0) the sensor representations
%   S.blocksize        - size of blocks used internally to split large
%                        continuous files [default ~20Mb]
%   S.updatehistory    - if 0 the history is not updated (use1ful for
%                        functions that use montage functionality.
%   S.prefix           - prefix for the output file (default - 'M')
%
% NOTE:
% montage are always defined based on the raw data on disk, i.e. discarding
% any curently applied online montage!
% Example: Data with 256 channels, online montage with a subset of 64
% channels. The montage must be based on the original 256 channels, not the
% "online" 64 ones.
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
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Robert Oostenveld, Stefan Kiebel, Christophe Phillips
% $Id: spm_eeg_montage.m 7169 2017-09-19 10:42:27Z vladimir $

SVNrev = '$Rev: 7169 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG montage'); spm('Pointer','Watch');

if ~isfield(S, 'montage') || isempty(S.montage)
    S.mode = 'clear';
    S.montage = 0;
elseif ~isfield(S, 'mode')
    S.mode      = 'write';
end
if ~isfield(S, 'blocksize'),     S.blocksize     = 655360;  end   %20 Mb
if ~isfield(S, 'keepothers'),    S.keepothers    = 1;       end
if ~isfield(S, 'keepsensors'),   S.keepsensors   = 1;       end
if ~isfield(S, 'updatehistory'), S.updatehistory = 1;       end
if ~isfield(S, 'prefix'),        S.prefix = 'M';            end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

if ~isequal(S.mode, 'write')
    D = D.montage('switch', 0);
end

%-Get montage
%--------------------------------------------------------------------------
if ischar(S.montage)
    montage = strmatch(S.montage, D.montage('getname', ':'));
    if isempty(montage)
        montage = load(S.montage);
    else
        error('Could not read montage file %s.', S.montage);
    end
    if ismember('montage', fieldnames(montage))
        montage = montage.montage;
    else
        error('Invalid montage file %s.', S.montage);
    end
else
    montage = S.montage;
end

if isequal(S.mode, 'write') && isnumeric(montage)
    if ~montage
        warning('The settings require applying the current montage. Nothing done');
        return;
    end
    montage = D.montage('getmontage', montage);
end

if ~isnumeric(montage)
    if ~all(isfield(montage, {'labelnew', 'labelorg', 'tra'})) || ...
            any([isempty(montage.labelnew), isempty(montage.labelorg), isempty(montage.tra)]) || ...
            length(montage.labelnew) ~= size(montage.tra, 1) || length(montage.labelorg) ~= size(montage.tra, 2)
        error('Insufficient or illegal montage specification.');
    end
    
    % There is inconsistent naming in FT functions.
    % montage.labelold = montage.labelorg; 
    
    % select and discard the columns that are empty if the corresponding
    % channels are not present in the data
    selempty   = find(all(montage.tra == 0, 1));
    
    [dum, sel] = setdiff(montage.labelorg(selempty), D.chanlabels);
    selempty   = selempty(sel);
    
    montage.tra(:, selempty)   = [];
    montage.labelorg(selempty) = [];
    
    if isfield(montage, 'chantypeold')
        montage.chantypeold(selempty) = [];
    end
    if isfield(montage, 'chanunitold')
        montage.chanunitold(selempty) = [];
    end
    
    % add columns for the channels that are not involved in the montage
    [add, ind] = setdiff(D.chanlabels, montage.labelorg);
    chlab      = D.chanlabels;
    ind        = sort(ind);
    add        = chlab(ind);
    
    m = size(montage.tra, 1);
    n = size(montage.tra, 2);
    k = length(add);
    if S.keepothers
        montage.tra((m+(1:k)), (n+(1:k))) = eye(k);
        montage.labelnew = cat(1, montage.labelnew(:), add(:));
        
        if isfield(montage, 'chantypenew')
            montage.chantypenew = lower(cat(1, montage.chantypenew(:), D.chantype(ind)'));
        end
        if isfield(montage, 'chanunitnew')
            montage.chanunitnew = cat(1, montage.chanunitnew(:), D.units(ind)');
        end
    else
        montage.tra = [montage.tra zeros(m, k)];
    end
    montage.labelorg = cat(1, montage.labelorg(:), add(:));
    
    if isfield(montage, 'chantypeold')
        montage.chantypeold = lower(cat(1, montage.chantypeold(:), D.chantype(ind)'));
    end
    if isfield(montage, 'chanunitold')
        montage.chanunitold = cat(1, montage.chanunitold(:), D.units(ind)');
    end    
    
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
    [selchan, selmont]  = spm_match_str(D.chanlabels, montage.labelorg);
    montage.tra         = montage.tra(:,selmont);
    montage.labelorg    = montage.labelorg(selmont);
    if isfield(montage, 'chantypeold')
        montage.chantypeold = montage.chantypeold(selmont);
    end
    if isfield(montage, 'chanunitold')
        montage.chanunitold = montage.chanunitold(selmont);
    end
end

isTF = strncmp(D.transformtype, 'TF', 2);

if isTF && ~isequal(S.mode, 'write')
    error('Online montages are not supported for TF data.');
end

switch S.mode
    case 'write'

        %-Generate new MEEG object with new filenames
        %------------------------------------------------------------------
        if ~isTF
            Dnew = clone(D, [S.prefix fname(D)], [m D.nsamples D.ntrials], 1);
            nblocksamples = floor(S.blocksize/max(D.nchannels, m));
        else
            Dnew = clone(D, [S.prefix fname(D)], [m D.nfrequencies D.nsamples D.ntrials], 1);
            nblocksamples = max(1, floor(S.blocksize/(D.nfrequencies*max(D.nchannels, m))));
        end
        
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
                if isTF
                    for f = 1:D.nfrequencies
                        temp = montage.tra * spm_squeeze(D(:, f, ((j-1)*nblocksamples +1) : (j*nblocksamples), i), [2, 4]);
                        Dnew(:, f, ((j-1)*nblocksamples +1) : (j*nblocksamples), i) = permute(shiftdim(temp, -1), [2, 1, 3]);
                    end
                else
                    Dnew(:, ((j-1)*nblocksamples +1) : (j*nblocksamples), i) = ...
                        montage.tra * spm_squeeze(D(:, ((j-1)*nblocksamples +1) : (j*nblocksamples), i), 3);
                end
                
                if D.ntrials == 1 && ismember(j, Ibar)
                    spm_progress_bar('Set', j);
                end
            end
            
            if mod(D.nsamples, nblocksamples) ~= 0
                if isTF
                    for f = 1:D.nfrequencies
                        temp = montage.tra * spm_squeeze(D(:, f, (nblocks*nblocksamples +1) : D.nsamples, i), [2, 4]);
                        Dnew(:, f, (nblocks*nblocksamples +1) : D.nsamples, i) = permute(shiftdim(temp, -1), [2, 1, 3]);
                    end
                else
                    Dnew(:, (nblocks*nblocksamples +1) : D.nsamples, i) = ...
                        montage.tra * squeeze(D(:, (nblocks*nblocksamples +1) : D.nsamples, i));
                end
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
        
        % If all the original channels contributing to a new channel have
        % the same units, transfer them to the new channel. This might be
        % wrong if the montage itself changes the units by scaling the data.
        if isfield(montage, 'chanunitnew')
            unitlist = montage.chanunitnew;
        else
            unitlist = repmat({'unknown'}, m, 1);
        end
        
        for i = 1:m
            if isequal(unitlist{i}, 'unknown')
                unit = unique(units(D, D.indchannel(montage.labelorg(~~abs(montage.tra(i, :))))));
                if numel(unit)==1
                    unitlist(i) = unit;
                end
            end
        end
        Dnew = units(Dnew, ':', unitlist);
        
        % Set channel types to default
        S1 = [];
        S1.task = 'defaulttype';
        S1.D = Dnew;
        S1.updatehistory = 0;
        Dnew = spm_eeg_prep(S1);
        
        if isfield(montage, 'chantypenew')
            ctype = Dnew.chantype;                        
            
            % Montage chatype overrides the default if specified
            for i = 1:m
                if ~isequal(montage.chantypenew{i}, 'unknown')
                    ctype(i) = upper(montage.chantypenew(i));
                end
            end
            
            Dnew = chantype(Dnew, ':', ctype);
        end
        
        %-Apply montage to sensors
        %------------------------------------------------------------------
        sensortypes = {'MEG', 'EEG'};
        for i = 1:length(sensortypes)
            if S.keepsensors
                sens = D.sensors(sensortypes{i});
                if ~isempty(sens) && ~isempty(intersect(sens.label, montage.labelorg))
                    sensmontage = montage;
                    [sel1, sel2] = spm_match_str(sens.label, sensmontage.labelorg);
                    sensmontage.labelorg = sensmontage.labelorg(sel2);
                    sensmontage.tra = sensmontage.tra(:, sel2);
                    selempty  = find(all(sensmontage.tra == 0, 2));
                    sensmontage.tra(selempty, :) = [];
                    sensmontage.labelnew(selempty) = [];
                    
                    if isfield(sensmontage, 'chantypeold')
                        sensmontage.chantypeold = sensmontage.chantypeold(sel2);
                    end
                    if isfield(sensmontage, 'chanunitold')
                        sensmontage.chanunitold = sensmontage.chanunitold(sel2);
                        
                        for c = 1:length(sensmontage.chanunitold)
                            if isequal(sensortypes{i}, 'MEG')
                                sensmontage.chanunitold{c} = strrep(sensmontage.chanunitold{c}, 'fT', 'T');
                            else
                                sensmontage.chanunitold{c} = strrep(sensmontage.chanunitold{c}, 'uV', 'V');
                            end
                            sensmontage.chanunitold{c} = strrep(sensmontage.chanunitold{c}, 'mm', 'm');
                        end
                        
                    end
                    if isfield(sensmontage, 'chantypenew')
                        sensmontage.chantypenew(selempty) = [];
                    end
                    if isfield(sensmontage, 'chanunitnew')
                        sensmontage.chanunitnew(selempty) = [];
                        
                        for c = 1:length(sensmontage.chanunitnew)
                            if isequal(sensortypes{i}, 'MEG')
                                sensmontage.chanunitnew{c} = strrep(sensmontage.chanunitnew{c}, 'fT', 'T');
                            else
                                sensmontage.chanunitnew{c} = strrep(sensmontage.chanunitnew{c}, 'uV', 'V');
                            end
                            sensmontage.chanunitnew{c} = strrep(sensmontage.chanunitnew{c}, 'mm', 'm');
                        end
                    end
                    
                    % Just remove known non-scalp channels to be on the
                    % safe side. 'Other' channels are not removed as they
                    % can be some kind of spatial components.
                    lblaux    = Dnew.chanlabels(Dnew.indchantype({'EOG', 'ECG', 'EMG', 'LFP', 'PHYS', 'ILAM', 'SRC'}));
                    [sel3, sel4] = spm_match_str(lblaux, sensmontage.labelnew);
                    sensmontage.tra(sel4, :) = [];
                    sensmontage.labelnew(sel4) = [];
                    
                    if S.keepothers
                        keepunused = 'yes';
                    else
                        keepunused = 'no';
                    end
                    
                    chanunitorig = sens.chanunit;
                    labelorg     = sens.label;
                    
                    sens = ft_apply_montage(sens, sensmontage, 'keepunused', keepunused, 'warning', false);
                    
                    if isequal(sensortypes{i}, 'MEG')
                        if isfield(sens, 'balance') && ~isequal(sens.balance.current, 'none')
                            balance = ft_apply_montage(getfield(sens.balance, sens.balance.current), sensmontage, 'keepunused', keepunused);
                        else
                            balance = sensmontage;
                        end
                        
                        sens.balance.custom = balance;
                        sens.balance.current = 'custom';
                    end
                    
                    
                    if isfield(sens, 'chanunit')
                        chanunit = sens.chanunit;
                    else
                        chanunit = repmat({'unknown'}, numel(sens.label), 1);
                    end
                    
                    % If all the original channels contributing to a new channel have
                    % the same units, transfer them to the new channel. This might be
                    % wrong if the montage itself changes the units by scaling the data.
                    for j = 1:numel(chanunit)
                        if isequal(chanunit{j}, 'unknown') %do not override units specified by montage
                            k = strmatch(sens.label{j}, sensmontage.labelnew, 'exact');
                            if ~isempty(k)
                                unit = unique(chanunitorig(sel1(~~abs(sensmontage.tra(k, :)))));
                                if numel(unit)==1
                                    chanunit(j) = unit;
                                end
                            else %channel was not in the montage, but just copied
                                k = strmatch(sens.label{j}, labelorg, 'exact');
                                chanunit(j) = chanunitorig(k);
                            end
                        end
                    end
                    
                    sens.chanunit = chanunit;
                end
                
                if ~isempty(sens) && ~isempty(intersect(sens.label, Dnew.chanlabels))
                    Dnew = sensors(Dnew, sensortypes{i}, sens);
                else
                    Dnew = sensors(Dnew, sensortypes{i}, []);
                end
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
        %------------------------------------------------------------------
        if ~isempty(Dnew.indchantype('EEG')) && (isempty(Dnew.sensors('EEG')) ||...
                (all(ismember({'spmnas', 'spmlpa', 'spmrpa'}, Dnew.fiducials.fid.label)) && ...
                ~isempty(setdiff(Dnew.chanlabels(Dnew.indchantype('EEG')), getfield(Dnew.sensors('EEG'), 'label')))))
            S1 = [];
            S1.task = 'defaulteegsens';
            S1.updatehistory = 0;
            S1.D = Dnew;
            
            Dnew = spm_eeg_prep(S1);
        end
        
        %-Create 2D positions for EEG by projecting the 3D positions to 2D
        %------------------------------------------------------------------
        if ~isempty(Dnew.indchantype('EEG')) && ~isempty(Dnew.sensors('EEG'))
            S1 = [];
            S1.task = 'project3D';
            S1.modality = 'EEG';
            S1.updatehistory = 0;
            S1.D = Dnew;
            
            Dnew = spm_eeg_prep(S1);
        end
        
        %-Create 2D positions for MEG  by projecting the 3D positions to 2D
        %------------------------------------------------------------------
        if ~isempty(Dnew.indchantype('MEG')) && ~isempty(Dnew.sensors('MEG'))
            S1 = [];
            S1.task = 'project3D';
            S1.modality = 'MEG';
            S1.updatehistory = 0;
            S1.D = Dnew;
            
            Dnew = spm_eeg_prep(S1);
        end
        
        %-Transfer the properties of channels not involved in the montage
        %------------------------------------------------------------------
        if ~isempty(add) && S.keepothers
            old_add_ind = D.indchannel(add);
            new_add_ind = Dnew.indchannel(add);
            
            Dnew = badchannels(Dnew, new_add_ind, badchannels(D, old_add_ind));
            Dnew = chantype(Dnew, new_add_ind, chantype(D, old_add_ind));
        end
    case 'switch'
        if isnumeric(montage)
            Dnew = D.montage('switch', montage);
        else
            Dnew = D.montage('add', montage);
        end
    case 'add'
        Dnew = D.montage('add', montage);
    case 'clear'
        Dnew = D.montage('clear');
end

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = check(Dnew);

if ~isfield(S, 'updatehistory') || S.updatehistory
    D = D.history('spm_eeg_montage', S);
end

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG montage: done'); spm('Pointer','Arrow');

% PROGRAMMER'S NOTES by CP
% Observed a weird behaviour of the montage method for the @meeg object
% under WinXP, matlab R2007b.
% 'montage' is initialized as a variable at "compilation" and then all
% calls like Nmont = montage(D,'getnumber'); are crashing.
% I therefore used the not so good looking Nmont = D.montage('getnumber');
% call. Same goes for the other calls to 'montage', even wih an extra
% argument.

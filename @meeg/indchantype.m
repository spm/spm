function ind = indchantype(this, types, flag)
% Method for getting channel indices based on labels and/or types
% FORMAT  ind = indchantype(this, types)
% this       - MEEG object
% channels   - string or cell array of strings may include
%             ('ALL', 'EEG', 'MEG', 'ECG', 'EOG' etc.)
% flag       - 'GOOD' or 'BAD' to include only good or bad channels
%              respectively (all are selected by default)
%              
% ind        - vector of channel indices matching labels
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: indchantype.m 6446 2015-05-22 14:10:56Z vladimir $

if ischar(types)    
    types = {types};
end

types = upper(types);
types = types(:)';

if ismember('ALL', types)
    ind = 1:nchannels(this);
else
    if ismember('FILTERED', types)
        types = [types, 'MEEG', 'REF', 'EOG', 'ECG', 'EMG', 'LFP', 'PHYS', 'ILAM', 'SRC'];
        types = setdiff(types, 'MEGCOMB');
    end
    
    if ismember('EOG', types)
        types = [types, 'VEOG', 'HEOG'];
    end
    
    if ismember('ECG', types)
        types = [types, 'EKG'];
    end
    
    if ismember('REF', types)
        types = [types, 'REFMAG', 'REFGRAD', 'REFPLANAR'];
    end
    
    if ismember('MEG', types)
        types = [types, 'MEGMAG', 'MEGGRAD'];
    end
    
    if ismember('MEGANY', types)
        types = [types, 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR'];
    end
    
    if ismember('MEEG', types)
        types = [types, 'EEG', 'MEG', 'MEGMAG', 'MEGCOMB', 'MEGGRAD', 'MEGPLANAR'];
    end
    
    ind = find(ismember(upper(chantype(this)), types));
end

if nargin > 2
    if strcmpi(flag, 'GOOD')
        ind = setdiff(ind, badchannels(this));
    elseif strcmpi(flag, 'BAD')
        ind = intersect(ind, badchannels(this));
    end
end

ind = sort(unique(ind));

ind = ind(:)'; % must be row to allow to use it as loop indices
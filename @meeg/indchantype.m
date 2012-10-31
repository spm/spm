function ind = indchantype(this, types)
% Method for getting channel indices based on labels and/or types
% FORMAT  ind = indchantype(this, types)
% this       - MEEG object
% channels   - string or cell array of strings may include
%             ('ALL', 'EEG', 'MEG', 'ECG', 'EOG' etc.)
%              If 'GOOD' or 'BAD' are included only good or bad channels
%              are selected accordingly.
%              
% ind        - vector of channel indices matching labels
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: indchantype.m 5025 2012-10-31 14:44:13Z vladimir $

if ischar(types)    
    types = {types};
end

types = upper(types);
types = types(:)';

if isequal(types, {'GOOD'})
    types = {'ALL', 'GOOD'};
end

if isequal(types, {'BAD'})
    types = {'ALL', 'BAD'};
end

if ismember('ALL', types)
    ind = 1:nchannels(this);
else
    if ismember('FILTERED', types)
        types = [types, 'MEEG', 'REF', 'EOG', 'ECG', 'EMG'];
    end
    
    if ismember('EOG', types)
        types = [types, 'VEOG', 'HEOG'];
    end
    
    if ismember('ECG', types)
        types = [types, 'EKG'];
    end
    
    if ismember('REF', types)
        types = [types, 'REFMAG', 'REFGRAD'];
    end
    
    if ismember('MEG', types)
        types = [types, 'MEGMAG', 'MEGGRAD'];
    end
    
    if ismember('MEEG', types)
        types = [types, 'EEG', 'MEG', 'MEGMAG', 'MEGGRAD', 'MEGPLANAR'];
    end
    
    ind = find(ismember(upper(chantype(this)), types));
end

if ismember('GOOD', types) && ~ismember('BAD', types)
        ind = setdiff(ind, badchannels(this));
end

if ismember('BAD', types) && ~ismember('GOOD', types)
        ind = intersect(ind, badchannels(this));
end
    
ind = sort(unique(ind));

ind = ind(:)'; % must be row to allow to use it as loop indices
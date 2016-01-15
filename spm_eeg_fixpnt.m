function data = spm_eeg_fixpnt(data, recurse)
% Helper function to replace pos by pnt
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_fixpnt.m 6683 2016-01-15 16:15:32Z guillaume $

if nargin==1
    recurse = 1;
end

if ~isa(data, 'struct')
    return;
end

if numel(data)>1
    % loop over all individual elements
    clear tmp
    for i=1:numel(data)
        % this is to prevent an "Subscripted assignment between dissimilar structures" error
        tmp(i) = spm_eeg_fixpnt(data(i));
    end
    data = tmp;
    clear tmp
    return
end

% replace pos by pnt
if isfield(data, 'pos')
    data.pnt = data.pos;
    data = rmfield(data, 'pos');
end

if recurse<3
    % recurse into substructures, not too deep
    fn = fieldnames(data);
    fn = setdiff(fn, {'cfg'}); % don't recurse into the cfg structure
    for i=1:length(fn)
        if isstruct(data.(fn{i}))
            data.(fn{i}) = spm_eeg_fixpnt(data.(fn{i}), recurse+1);
        end
    end
end

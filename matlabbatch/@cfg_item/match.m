function sts = match(item, spec)

% MATCH function for cfg trees
% This function is an implementation of find to search the cfg tree for
% certain entries.
%
% sts = match(item, spec)
% Spec must be a cell array of struct arrays with one or more fields. Each
% struct must contain two fields - 'name' and 'value'.
% An item matches, if it has a field with the specified field name and the
% contents of this field equals the contents of spec.value. If the field
% name is 'class', an item matches, if its class name is equal to
% spec.value.
% Matches within each struct array are OR-concatenated, while matches
% between struct arrays are AND-concatenated.
% An empty spec always matches.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: match.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

% match an empty spec
sts = true;

for k = 1:numel(spec)
    % Assume no match
    sts = false;
    for l = 1:numel(spec{k})
        if strcmp(spec{k}(l).name, 'class')
            sts = strcmp(spec{k}(l).value, class(item));
        else
            try
                % This may fail if there is no field with this name in item
                sts = isequal(spec{k}(l).value, ...
                    subsref(item, substruct('.', spec{k}(l).name)));
            end;
        end;
        if sts
            % OR: success on first match
            break;
        end;
    end;
    if ~sts
        % AND: fail on first non-match
        break;
    end;
end;

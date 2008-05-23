function sts = match(item, spec)

% function sts = match(item, spec)
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
% Special matching rules for cfg_files apply to the .filter field, if
% both item.filter and spec.value are one of the special types 'any',
% 'image', 'nifti', 'mat', 'xml', 'batch', 'dir':
% A .filter 'any' matches any spec.value. All other filters only match if
% strcmpi(item.filter,spec.value) is true. Currently, 'nifti' and 'image'
% filters are treated as equivalent.
% Checking the equivalence of two regular expressions is a demanding
% task. Therefore, no matching is performed if item.filter or spec.value
% are regular expressions and this match will always be true.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: match.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

% match an empty spec
sts = true;

specflt = {'image','nifti','mat','mattxt','batch','xml','any','dir'};
for k = 1:numel(spec)
    % Assume no match
    sts = false;
    for l = 1:numel(spec{k})
        switch spec{k}(l).name,
            % don't try any matching of regexp filters
            case 'filter',
                if strcmpi(item.filter,'nifti')
                    ifilter = 'image';
                else
                    ifilter = item.filter;
                end;
                if strcmpi(spec{k}(l).value,'nifti')
                    sfilter = 'image';
                else
                    sfilter = spec{k}(l).value;
                end;                
                if strcmpi(ifilter,'any')
                    sts = true;
                elseif any(strcmpi(sfilter,specflt)) && any(strcmpi(ifilter,specflt))
                    sts = strcmpi(sfilter,ifilter);
                else
                    sts = true;
                end;
            case 'ufilter',
                sts = true;
            case 'num',
                sts = true;
            case 'class'
                sts = strcmpi(spec{k}(l).value,class(item));
            otherwise
                spec1{1}(1) = spec{k}(l);
                sts = match(item.cfg_item, spec1);
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

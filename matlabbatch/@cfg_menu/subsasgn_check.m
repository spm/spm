function [sts, val] = subsasgn_check(item,subs,val)

% function [sts, val] = subsasgn_check(item,subs,val)
% Perform assignment checks for .val, .labels and .values field. 
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

% Checks could include a check whether the number of elements in labels
% and values match, but this would require a method to add a label and a
% value simultaneously. 
% Also, a check whether .val is indeed in .values could be performed, if
% the .val field is filled after .values (which is not the case for
% current auto-generated code).

sts = true;
switch subs(1).subs
    case {'val'}
        sts = iscell(val) && (isempty(val) || numel(val) == 1);
        if ~sts
            cfg_message('matlabbatch:checkval', ...
                    '%s: Value must be a cell with zero or one elements.', subsasgn_checkstr(item,subs));
        end        
    case {'labels'}
        sts = iscell(val) && (isempty(val) || iscellstr(val));
        if ~sts
            cfg_message('matlabbatch:check:labels', ...
                    '%s: Value must be a cell array of strings.', subsasgn_checkstr(item,subs));
        end
    case {'values'}
        sts = iscell(val);
        if ~sts
            cfg_message('matlabbatch:check:values', ...
                    '%s: Value must be a cell array of values.', subsasgn_checkstr(item,subs));
        end
end

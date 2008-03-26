function [sts val] = subsasgn_check(item,subs,val)

% function [sts val] = subsasgn_check(item,subs,val)
% Perform assignment checks for .labels and .values field. Checks could
% include a check whether the number of elements in labels and values
% match, but this would require a method to add a label and a value
% simultaneously. 
% Also, a check whether .val is indeed in .values could be performed.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1245 2008-03-26 10:32:25Z volkmar $

rev = '$Rev: 1245 $';

sts = true;
checkstr = sprintf('Item ''%s'', field ''%s''', subsref(item,substruct('.','name')), subs(1).subs);
switch subs(1).subs
    case {'labels'}
	sts = iscell(val) && (isempty(val) || iscellstr(val));
        if ~sts
            warning('matlabbatch:cfg_menu:subsasgn_check:labels', ...
                    '%s: Value must be a cell array of strings.', checkstr);
        end;
    case {'values'}
        sts = iscell(val);
        if ~sts
            warning('matlabbatch:cfg_menu:subsasgn_check:values', ...
                    '%s: Value must be a cell array of values.', checkstr);
        end;
end;

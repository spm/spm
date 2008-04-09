function [sts val] = subsasgn_check(item,subs,val)

% function [sts val] = subsasgn_check(item,subs,val)
% Check whether .prog, .vout and .vfiles are functions or function
% handles and whether dependencies are cfg_dep objects.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1322 2008-04-09 08:37:52Z volkmar $

rev = '$Rev: 1322 $';

sts = true;
checkstr = sprintf('Item ''%s'', field ''%s''', subsref(item,substruct('.','name')), subs(1).subs);
switch subs(1).subs
    case {'prog', 'vout', 'vfiles'}
        sts = subsasgn_check_funhandle(val);
        if ~sts
            warning('matlabbatch:cfg_exbranch:subsasgn_check', ...
                    ['%s: Value must be a function or function handle on ' ...
                     'MATLAB path.'], checkstr);
        end;
    case {'sdeps', 'tdeps', 'sout'}
        sts = isempty(val) || isa(val, 'cfg_dep');
        if ~sts
            warning('matlabbatch:cfg_exbranch:subsasgn_check', ...
                    '%s: Value must be a cfg_dep dependency object.', checkstr);
        end;
end;

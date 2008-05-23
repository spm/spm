function [sts, val] = subsasgn_check(item,subs,val)

% function [sts, val] = subsasgn_check(item,subs,val)
% Check whether .prog, .vout and .vfiles are functions or function
% handles and whether dependencies are cfg_dep objects.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

sts = true;
switch subs(1).subs
    case {'prog', 'vout', 'vfiles'}
        sts = subsasgn_check_funhandle(val);
        if ~sts
            warning('matlabbatch:cfg_exbranch:subsasgn_check', ...
                    ['%s: Value must be a function or function handle on ' ...
                     'MATLAB path.'], subsasgn_checkstr(item,subs));
        end;
    case {'sdeps', 'tdeps', 'sout'}
        sts = isempty(val) || isa(val, 'cfg_dep');
        if ~sts
            warning('matlabbatch:cfg_exbranch:subsasgn_check', ...
                    '%s: Value must be a cfg_dep dependency object.', subsasgn_checkstr(item,subs));
        end;
end;

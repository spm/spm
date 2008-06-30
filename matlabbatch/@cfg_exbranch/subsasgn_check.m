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
% $Id: subsasgn_check.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

sts = true;
switch subs(1).subs
    case {'prog', 'vout', 'vfiles'}
        sts = subsasgn_check_funhandle(val);
        if ~sts
            cfg_message('matlabbatch:check:funhandle', ...
                    ['%s: Value must be a function or function handle on ' ...
                     'MATLAB path.'], subsasgn_checkstr(item,subs));
        end
    case {'sdeps', 'tdeps', 'sout'}
        sts = isempty(val) || isa(val, 'cfg_dep');
        if ~sts
            cfg_message('matlabbatch:check:dep', ...
                    '%s: Value must be a cfg_dep dependency object.', subsasgn_checkstr(item,subs));
        end
end

function [sts val] = subsasgn_check(item,subs,val)

% function [sts val] = subsasgn_check(item,subs,val)
% Check whether .prog, .vout and .vfiles are functions or function handles.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

sts = true;

[sts val] = subsasgn_check(item.cfg_branch,subs,val);
if ~sts
    return;
end;
% check, whether arguments for 'val' are cfg_items
switch subs(1).subs
    case {'prog', 'vout', 'vfiles'}
        sts = subsasgn_check_funhandle(val);
end;

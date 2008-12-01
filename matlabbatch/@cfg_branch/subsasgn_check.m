function [sts, val] = subsasgn_check(item,subs,val)

% function [sts, val] = subsasgn_check(item,subs,val)
% Check whether type of val conforms to configuration tree specification.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 2512 2008-12-01 13:21:29Z volkmar $

rev = '$Rev: 2512 $'; %#ok

sts = true;

% check, whether arguments for 'val' are cfg_items
switch subs(1).subs
    case {'val'}
        sts = ~cfg_get_defaults('cfg_item.checkval') || subsasgn_check_valcfg(subs,val,[0 Inf]);
end;

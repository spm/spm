function [sts, val] = subsasgn_check(item,subs,val)

% function [sts, val] = subsasgn_check(item,subs,val)
% Check assignments to item.values and item.val. 
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

sts = true;
switch subs(1).subs
    case {'values'}
	sts = subsasgn_check_valcfg(subs,val,[0 Inf]);
    case {'val'}
	% Check maximum number of elements
	sts = subsasgn_check_valcfg(subs,val,[0 1]);
	% Could also check whether added element is one from 'values' list
end;

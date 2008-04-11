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
% $Id: subsasgn_check.m 1373 2008-04-11 14:24:03Z spm $

rev = '$Rev: 1373 $';

sts = true;

% check, whether arguments for 'val' are cfg_items
switch subs(1).subs
    case {'val'}
    sts = subsasgn_check_valcfg(subs,val,[0 Inf]);
end;

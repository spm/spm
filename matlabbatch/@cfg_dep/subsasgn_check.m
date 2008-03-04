function [sts val] = subsasgn_check(item,subs,val)

% function [sts val] = subsasgn_check(item,subs,val)
% Do a check for proper assignments of values to fields. This routine
% will be called for derived objects from @cfg_item/subsasgn with the
% original object as first argument and the proposed subs and val fields
% before an assignment is made. It is up to each derived class to
% implement assignment checks for both its own fields and fields
% inherited from cfg_item.
% 
% This routine is a fallback placeholder in cfg_item if a derived class
% does not need to implement its own checks and always returns true.
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
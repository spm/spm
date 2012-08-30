function fn = fieldnames(item)

% function fn = fieldnames(item)
% Return a list of all (inherited and non-inherited) field names.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: fieldnames.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

fn1 = fieldnames(item.cfg_item);
fn2 = mysubs_fields;

fn = [fn1(:); fn2(:)]';
function fn = fieldnames(item)

% Return a list of all (inherited and non-inherited) field names.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: fieldnames.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

fn = mysubs_fields;
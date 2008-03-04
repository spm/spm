function fnames = subs_fields(item)

% This function works as a "class-based switch" to return the value of
% the private mysubs_fields function for the appropriate class. It is
% identical for each class, but it has to be in the class directory to
% access the proper private function mysubs_fields.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subs_fields.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

fnames = mysubs_fields;
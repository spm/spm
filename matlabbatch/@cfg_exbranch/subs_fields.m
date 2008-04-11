function fnames = subs_fields(item)

% function fnames = subs_fields(item)
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
% $Id: subs_fields.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

fnames = mysubs_fields;
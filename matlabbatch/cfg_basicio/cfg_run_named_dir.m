function out = cfg_run_named_dir(job)

% Return selected dirs as separate output arguments.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_named_dir.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

out.dirs = job.dirs;

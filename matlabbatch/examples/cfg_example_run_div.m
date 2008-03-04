function out = cfg_example_run_div(job)
% Example function that returns the mod and rem of two numbers given in
% job.a and job.b in out.mod and out.rem.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_run_div.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

out.mod = mod(job.a, job.b);
out.rem = mod(job.a, job.b);
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
% $Id: cfg_example_run_div.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

out.mod = mod(job.a, job.b);
out.rem = mod(job.a, job.b);
function out = cfg_example_run_cumsum2(job)
% Example function that returns the cumulative sum of an vector given in
% job.a in out. The output is referenced as out(:), this is defined in
% cfg_example_vout_cumsum1.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_run_cumsum2.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

% The harvested cfg_repeat will return a cell array of numbers. These need
% to be cat'ed together to be useful as input to cumsum.
out.cs = cumsum(cat(1,job.a{:}));
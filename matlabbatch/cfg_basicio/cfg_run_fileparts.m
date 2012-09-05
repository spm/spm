function out = cfg_run_fileparts(job)

% Run fileparts on a list of files.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_fileparts.m 4899 2012-09-05 13:44:17Z volkmar $

rev = '$Rev: 4899 $'; %#ok

[out.p out.n out.e] = cellfun(@fileparts, job.files, 'UniformOutput', false);
out.up = unique(out.p);
function out = cfg_run_file_fplist(job)

% function out = cfg_run_file_fplist(job)
%
% Select files non-interactively using cfg_getfile('FPList',...)
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_file_fplist.m 1896 2008-07-09 08:21:36Z volkmar $

rev = '$Rev: 1896 $'; %#ok

[out.files out.dirs] = cfg_getfile('FPList', job.dir{1}, job.filter);

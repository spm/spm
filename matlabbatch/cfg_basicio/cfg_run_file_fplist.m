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
% $Id: cfg_run_file_fplist.m 1203 2008-03-12 09:29:49Z volkmar $

rev = '$Rev: 1203 $';

[files dirs] = cfg_getfile('FPList', job.dir{1}, job.filter);
out.files = cellstr(files);
out.dirs  = cellstr(dirs);
function out = cfg_run_file_fplist(job)

% function out = cfg_run_file_fplist(job)
%
% Select files non-interactively using cfg_getfile('FPList',...) or
% cfg_getfile('FPListRec',...).
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_file_fplist.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

[out.files, out.dirs] = cfg_getfile(job.rec, job.dir{1}, job.filter);

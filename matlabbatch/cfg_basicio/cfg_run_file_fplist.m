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
% $Id: cfg_run_file_fplist.m 6918 2016-11-02 14:33:11Z guillaume $

rev = '$Rev: 6918 $'; %#ok

[out.files, out.dirs] = cfg_getfile(job.rec, job.dir{1}, job.filter);
if numel(out.files) < 2, plural = ''; else plural = 's'; end
if strcmpi(job.rec,'FPListRec'), rec = ' (and subdirectories)'; else rec = ''; end
fprintf('Found %d file%s matching filter ''%s'' in %s%s.\n',...
    numel(out.files),plural,job.filter,job.dir{1},rec);

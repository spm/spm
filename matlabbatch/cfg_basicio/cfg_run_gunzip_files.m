function out = cfg_run_gunzip_files(job)
% Run gunzip on a list of files.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_gunzip_files.m 6628 2015-12-04 08:20:51Z volkmar $

rev = '$Rev: 6628 $'; %#ok

if isempty(job.outdir) || isempty(job.outdir{1})
    out = gunzip(job.files);
else
    out = gunzip(job.files, job.outdir{1});
end
if ~job.keep
    delete(job.files{:});
end
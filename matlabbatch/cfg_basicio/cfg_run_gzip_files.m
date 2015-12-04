function out = cfg_run_gzip_files(job)
% Run gzip on a list of files.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_gzip_files.m 6627 2015-12-04 08:20:49Z volkmar $

rev = '$Rev: 6627 $'; %#ok

if isempty(job.outdir)
    out = gzip(job.files);
else
    out = gzip(job.files, job.outdir{1});
end
if ~job.keep
    delete(job.files{:});
end
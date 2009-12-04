function out = spm_run_cat(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_run_cat.m 3613 2009-12-04 18:47:59Z guillaume $

V     = strvcat(job.vols{:});
dt    = job.dtype;
fname = job.name;

V4    = spm_file_merge(V,fname,dt);

out.mergedfile = {V4(1).fname};

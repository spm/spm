function out = spm_run_reorient(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_reorient.m 1185 2008-03-04 16:31:21Z volkmar $

job = varargin{1};
if isfield(job.transform,'transprm')
    job.transform.transM = spm_matrix(job.transform.transprm);
end;
spm_progress_bar('Init', numel(job.srcfiles), 'Reorient', 'Images completed');
for k = 1:numel(job.srcfiles)
    M = spm_get_space(job.srcfiles{k});
    spm_get_space(job.srcfiles{k},job.transform.transM*M);
    spm_progress_bar('Set',k);
end;
spm_progress_bar('Clear');
out.files = job.srcfiles;
return;

function out = convert_minc(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_minc.m 2312 2008-10-07 17:02:46Z volkmar $

out.files = cell(size(job.data));
for i=1:length(job.data),
    spm_mnc2nifti(job.data{i},job.opts);
    [pth,nam,ext,num] = spm_fileparts(job.data{i});
    out.files{i} = fullfile(pth,[nam job.opts.ext num]);
end;
return;

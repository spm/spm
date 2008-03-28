function out = spm_run_smooth(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_smooth.m 1274 2008-03-28 16:22:43Z volkmar $

job = varargin{1};
out.files = cell(size(job.data));
spm_progress_bar('Init',numel(job.data),'Smoothing','Volumes Complete');
for i = 1:numel(job.data)
        [pth,nam,ext,num] = spm_fileparts(job.data{i});
        out.files{i} = fullfile(pth,[job.prefix nam ext num]);
        spm_smooth(job.data{i},out.files{i},job.fwhm,job.dtype);
        spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

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

% $Id: spm_run_smooth.m 1185 2008-03-04 16:31:21Z volkmar $

job = varargin{1};

P     = strvcat(job.data);
s     = job.fwhm;
dtype = job.dtype;
n     = size(P,1);

spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
for i = 1:n
        Q = deblank(P(i,:));
        [pth,nam,ext,num] = spm_fileparts(deblank(Q));
        out.files{i} = fullfile(pth,[job.prefix nam ext num]);
        spm_smooth(Q,U,s,dtype);
        spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

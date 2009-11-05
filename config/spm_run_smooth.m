function out = spm_run_smooth(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_smooth.m 3534 2009-11-05 12:34:21Z guillaume $

job       = varargin{1};
out.files = cell(size(job.data));
spm_progress_bar('Init',numel(job.data),'Smoothing','Volumes Complete');
for i = 1:numel(job.data)
    [pth,nam,ext,num] = spm_fileparts(job.data{i});
    out.files{i}      = fullfile(pth,[job.prefix nam ext num]);
    spm_smooth(job.data{i},out.files{i},job.fwhm,job.dtype);
    if job.im
        vi = spm_read_vols(spm_vol(job.data{i}),1);
        v  = spm_vol(out.files{i});
        vo = spm_read_vols(v);
        if spm_type(v.dt(1),'nanrep')
            vo(isnan(vi)) = NaN;
        else
            vo(isnan(vi)) = 0;
        end
        spm_write_vol(v,vo);
    end
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

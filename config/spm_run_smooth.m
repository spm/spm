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

% $Id: spm_run_smooth.m 3915 2010-06-02 17:09:10Z guillaume $

job       = varargin{1};
out.files = cell(size(job.data));
spm_progress_bar('Init',numel(job.data),'Smoothing','Volumes Complete');
for i = 1:numel(job.data)
    [pth,nam,ext,num] = spm_fileparts(job.data{i});
    out.files{i}      = fullfile(pth,[job.prefix nam ext num]);
    spm_smooth(job.data{i},out.files{i},job.fwhm,job.dtype);
    if job.im
        vi = spm_vol(job.data{i});
        vo = spm_vol(out.files{i});
        for j=1:numel(vi)
            vvi = spm_read_vols(vi(j),1);
            vvo = spm_read_vols(vo(j));
            if spm_type(vo(j).dt(1),'nanrep')
                vvo(isnan(vvi)) = NaN;
            else
                vvo(isnan(vvi)) = 0;
            end
            spm_write_vol(vo(j),vvo);
        end
    end
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

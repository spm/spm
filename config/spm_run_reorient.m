function out = spm_run_reorient(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2006-2012 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_run_reorient.m 4649 2012-02-06 15:55:04Z guillaume $

job = varargin{1};
if isfield(job.transform,'transprm')
    job.transform.transM = spm_matrix(job.transform.transprm);
end
spm_progress_bar('Init', numel(job.srcfiles), 'Reorient', 'Images completed');
if isempty(job.prefix)
    for k = 1:numel(job.srcfiles)
        M = spm_get_space(job.srcfiles{k});
        spm_get_space(job.srcfiles{k},job.transform.transM*M);
        spm_progress_bar('Set',k);
    end
    out.files = job.srcfiles;
else
    out.files = cell(size(job.srcfiles));
    for k = 1:numel(job.srcfiles)
        V       = spm_vol(job.srcfiles{k});
        X       = spm_read_vols(V);
        V.mat   = job.transform.transM * V.mat;
        V.fname = spm_file(V.fname, 'prefix',job.prefix);
        spm_write_vol(V,X);
        out.files{k} = V.fname;
        spm_progress_bar('Set',k);
    end
end
spm_progress_bar('Clear');

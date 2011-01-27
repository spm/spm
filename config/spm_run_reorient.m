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

% $Id: spm_run_reorient.m 4177 2011-01-27 14:35:32Z volkmar $

job = varargin{1};
if isfield(job.transform,'transprm')
    job.transform.transM = spm_matrix(job.transform.transprm);
end;
spm_progress_bar('Init', numel(job.srcfiles), 'Reorient', 'Images completed');
if isempty(job.prefix)
    for k = 1:numel(job.srcfiles)
        M = spm_get_space(job.srcfiles{k});
        spm_get_space(job.srcfiles{k},job.transform.transM*M);
        spm_progress_bar('Set',k);
    end;
    out.files = job.srcfiles;
else
    out.files = cell(size(job.srcfiles));
    for k = 1:numel(job.srcfiles)
        V = spm_vol(job.srcfiles{k});
        X = spm_read_vols(V);
        [p n e v] = spm_fileparts(V.fname);
        V.mat = job.transform.transM*V.mat;
        V.fname = fullfile(p, [job.prefix n e v]);
        spm_write_vol(V,X);
        out.files{k} = V.fname;
        spm_progress_bar('Set',k);
    end
end
spm_progress_bar('Clear');
return;

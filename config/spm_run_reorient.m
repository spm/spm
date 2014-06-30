function out = spm_run_reorient(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2006-2014 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_run_reorient.m 6078 2014-06-30 18:10:33Z guillaume $

job = varargin{1};
if isfield(job.transform,'transprm')
    job.transform.transM = spm_matrix(job.transform.transprm);
elseif isfield(job.transform,'transF')
    load(char(job.transform.transF), 'M');
    job.transform.transM = M;
end

K = numel(job.srcfiles);
spm_progress_bar('Init', K, 'Reorient', 'Images completed');
if isempty(job.prefix)
    
    % read and write separately, so duplicates get harmlessly overwritten
    MM = zeros(4, 4, K);
    for k = 1:K
        MM(:, :, k) = spm_get_space(job.srcfiles{k});
    end
    for k = 1:K
        spm_get_space(job.srcfiles{k}, job.transform.transM * MM(:, :, k));
        spm_progress_bar('Set',k);
    end
    out.files = job.srcfiles;
    
else
    
    out.files = cell(size(job.srcfiles));
    for k = 1:K
        V       = spm_vol(job.srcfiles{k});
        X       = spm_read_vols(V);
        V.mat   = job.transform.transM * V.mat;
        V.fname = spm_file(V.fname, 'prefix',job.prefix);
        spm_write_vol(V,X);
        out.files{k} = [V.fname ',' num2str(V.n(1))];
        spm_progress_bar('Set',k);
    end
    
end
spm_progress_bar('Clear');

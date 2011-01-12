function out = spm_run_coreg_estimate(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_coreg_estimate.m 4157 2011-01-12 18:19:38Z guillaume $

job = varargin{1};

x  = spm_coreg(char(job.ref), char(job.source), job.eoptions);

M  = spm_matrix(x);
PO = [job.source(:); job.other(:)];
PO = PO(~cellfun('isempty',PO));
MM = zeros(4,4,numel(PO));
for j=1:numel(PO)
    MM(:,:,j) = spm_get_space(PO{j});
end
for j=1:numel(PO)
    spm_get_space(PO{j}, M\MM(:,:,j));
end

out.cfiles = PO;

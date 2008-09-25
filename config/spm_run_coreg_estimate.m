function out = spm_run_coreg_estimate(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_coreg_estimate.m 2187 2008-09-25 11:50:50Z volkmar $

job = varargin{1};
%disp(job);
%disp(job.eoptions);

x  = spm_coreg(strvcat(job.ref), strvcat(job.source),job.eoptions);
M  = inv(spm_matrix(x));
PO = strvcat(strvcat(job.source),strvcat(job.other));
MM = zeros(4,4,size(PO,1));
for j=1:size(PO,1),
    MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
end;
for j=1:size(PO,1),
    spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
end;
out.cfiles = cellstr(PO);
return;

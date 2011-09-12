function out = spm_run_coreg(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_coreg.m 4480 2011-09-12 14:34:50Z guillaume $


if ~isfield(job,'other') || isempty(job.other{1}), job.other = {}; end
PO = [job.source(:); job.other(:)];

%-Coregister
%--------------------------------------------------------------------------
if isfield(job,'eoptions')
    x  = spm_coreg(char(job.ref), char(job.source), job.eoptions);
    
    M  = spm_matrix(x);
    iM = inv(M);
    for j=1:numel(PO)
        MM = spm_get_space(PO{j});
        spm_get_space(PO{j}, iM*MM);
    end
end

%-Reslice
%--------------------------------------------------------------------------
if isfield(job,'roptions')
    P            = char(job.ref{:},job.source{:},job.other{:});
    flags.mask   = job.roptions.mask;
    flags.mean   = 0;
    flags.interp = job.roptions.interp;
    flags.which  = 1;
    flags.wrap   = job.roptions.wrap;
    flags.prefix = job.roptions.prefix;

    spm_reslice(P, flags);
end

%-Dependencies
%--------------------------------------------------------------------------
if isfield(job,'eoptions')
    out.cfiles   = PO;
    out.M        = M;
end
if isfield(job,'roptions')
    out.rfiles   = spm_file(PO, 'prefix',job.roptions.prefix);
end

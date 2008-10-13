function out = spm_run_normalise_estimate(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_normalise_estimate.m 2335 2008-10-13 18:49:25Z christophe $

job    = varargin{1};
o      = job.eoptions;
eflags = struct(...
    'smosrc', o.smosrc,...
    'smoref', o.smoref,...
    'regtype',o.regtype,...
    'cutoff', o.cutoff,...
    'nits',   o.nits,...
    'reg',    o.reg);

for i=1:numel(job.subj),
    [pth,nam,ext,ind] = spm_fileparts(strvcat(job.subj(i).source{:}));
    out(i).params{1}  = fullfile(pth,[nam,'_sn.mat']);
    spm_normalise(strvcat(job.eoptions.template{:}),...
        strvcat(job.subj(i).source{:}), out(i).params{1},...
        strvcat(job.eoptions.weight), strvcat(job.subj(i).wtsrc), eflags);
end;
return;

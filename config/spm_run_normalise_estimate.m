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

% $Id: spm_run_normalise_estimate.m 4152 2011-01-11 14:13:35Z volkmar $

job    = varargin{1};
eflags = rmfield(job.eoptions,{'template','weight'});

out = repmat(struct('params',{''}),size(job.subj));
for i=1:numel(job.subj),
    [pth,nam] = spm_fileparts(char(job.subj(i).source));
    out(i).params  = {fullfile(pth,[nam,'_sn.mat'])};
    spm_normalise(char(job.eoptions.template),...
        char(job.subj(i).source), out(i).params{1},...
        char(job.eoptions.weight), char(job.subj(i).wtsrc), eflags);
end;
return;

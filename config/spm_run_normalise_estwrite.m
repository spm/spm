function out = spm_run_normalise_estwrite(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_normalise_estwrite.m 4152 2011-01-11 14:13:35Z volkmar $

job    = varargin{1};
eflags = rmfield(job.eoptions,{'weight','template'});
rflags = job.roptions;

out = repmat(struct('params',{''},'files',{''}),size(job.subj));
for i=1:numel(job.subj),
    [ pth,nam ] = spm_fileparts(char(job.subj(i).source));
    out(i).params  = {fullfile(pth,[nam,'_sn.mat'])};
    spm_normalise(char(job.eoptions.template),...
        char(job.subj(i).source), out(i).params{1},...
        char(job.eoptions.weight), char(job.subj(i).wtsrc), eflags);
    spm_write_sn(char(job.subj(i).resample), out(i).params{1}, rflags);
    out(i).files = cell(size(job.subj(i).resample));
    for j=1:numel(job.subj(i).resample),
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).resample{j});
        out(i).files{j}   = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
    end;
end;
return;

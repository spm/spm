function out = spm_run_coreg_reslice(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_coreg_reslice.m 4152 2011-01-11 14:13:35Z volkmar $

job = varargin{1};

P            = char([job.ref(:);job.source(:)]);
flags.mask   = job.roptions.mask;
flags.mean   = 0;
flags.interp = job.roptions.interp;
flags.which  = 1;
flags.wrap   = job.roptions.wrap;
flags.prefix = job.roptions.prefix;

spm_reslice(P,flags);
out.rfiles = cell(size(job.source));
for i=1:numel(job.source),
    [pth,nam,ext,num] = spm_fileparts(job.source{i});
    out.rfiles{i} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
end;

return;

function out = spm_run_normalise_write(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_normalise_write.m 2312 2008-10-07 17:02:46Z volkmar $

job    = varargin{1};
o      = job.roptions;
rflags = struct(...
    'preserve',o.preserve,...
    'bb',      o.bb,...
    'vox',     o.vox,...
    'interp',  o.interp,...
    'wrap',    o.wrap,...
        'prefix',  o.prefix);

for i=1:numel(job.subj),
    spm_write_sn(strvcat(job.subj(i).resample{:}),...
        strvcat(job.subj(i).matname{:}),rflags);
    res = job.subj(i).resample;
    out(i).files = cell(size(res));
    for j=1:numel(res),
        [pth,nam,ext,num] = spm_fileparts(res{j});
        out(i).files{j} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
    end;
end;
return;

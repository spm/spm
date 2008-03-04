function out = spm_run_realign_reslice(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_realign_reslice.m 1185 2008-03-04 16:31:21Z volkmar $

job          = varargin{1};
P            = strvcat(job.data);
flags.mask   = job.roptions.mask;
flags.mean   = job.roptions.which(2);
flags.interp = job.roptions.interp;
flags.which  = job.roptions.which(1);
flags.wrap   = job.roptions.wrap;
flags.prefix = job.roptions.prefix;
spm_reslice(P,flags);
switch job.roptions.which(1),
case 1,
    out.rfiles = cell(numel(job.data)-1,1);
    for i=1:length(out.rfiles),
        [pth,nam,ext,num] = spm_fileparts(job.data{i+1});
        out.rfiles{i} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
    end;
otherwise,
    out.rfiles = cell(numel(job.data),1);
    for i=1:length(out.rfiles),
        [pth,nam,ext,num] = spm_fileparts(job.data{i});
        out.rfiles{i} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
    end;
end;
if job.roptions.which(2),
    [pth,nam,ext,num] = spm_fileparts(job.data{1});
    out.rmean{1} = fullfile(pth,['mean', nam, ext, num]);
end;
return;

function out = spm_run_coreg_estwrite(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_coreg_estwrite.m 1185 2008-03-04 16:31:21Z volkmar $

job = varargin{1};
%disp(job);
%disp(job.eoptions);
%disp(job.roptions);

job.ref    = strvcat(job.ref);
job.source = strvcat(job.source);
job.other  = strvcat(job.other);

x  = spm_coreg(job.ref, job.source,job.eoptions);
M  = inv(spm_matrix(x));
PO = strvcat(job.source,job.other);
MM = zeros(4,4,size(PO,1));
for j=1:size(PO,1),
        MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
end;
for j=1:size(PO,1),
        spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
end;

P            = strvcat(job.ref,job.source,job.other);
flags.mask   = job.roptions.mask;
flags.mean   = 0;
flags.interp = job.roptions.interp;
flags.which  = 1;
flags.wrap   = job.roptions.wrap;
flags.prefix = job.roptions.prefix;

spm_reslice(P,flags);

out.cfiles = cellstr(strvcat(job.source{:}, job.other{:}));
for i=1:numel(out.cfiles),
    [pth,nam,ext,num] = spm_fileparts(out.cfiles{i});
    out.rfiles{i} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
end;
return;

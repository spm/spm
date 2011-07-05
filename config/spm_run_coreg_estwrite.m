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

% $Id: spm_run_coreg_estwrite.m 4380 2011-07-05 11:27:12Z volkmar $

job = varargin{1};
if isempty(job.other{1})
    job.other = {};
end

x  = spm_coreg(char(job.ref), char(job.source),job.eoptions);

M  = spm_matrix(x);
PO = [job.source(:); job.other(:)];
MM = zeros(4,4,numel(PO));
for j=1:numel(PO),
    MM(:,:,j) = spm_get_space(PO{j});
end
for j=1:numel(PO),
    spm_get_space(PO{j}, M\MM(:,:,j));
end

P            = char(job.ref{:},job.source{:},job.other{:});
flags.mask   = job.roptions.mask;
flags.mean   = 0;
flags.interp = job.roptions.interp;
flags.which  = 1;
flags.wrap   = job.roptions.wrap;
flags.prefix = job.roptions.prefix;

spm_reslice(P,flags);

out.cfiles = PO;
out.M      = M;
out.rfiles = cell(size(out.cfiles));
for i=1:numel(out.cfiles),
    [pth,nam,ext,num] = spm_fileparts(out.cfiles{i});
    out.rfiles{i} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
end;
return;


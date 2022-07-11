function out = spm_run_st(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging

Seq = job.so;
TR  = job.tr;
TA  = job.ta;
nslices   = job.nslices;
refslice  = job.refslice;
timing(2) = TR - TA;
timing(1) = TA / (nslices -1);

P = cell(size(job.scans));
for i = 1:numel(job.scans)
    P{i} = char(job.scans{i});
end

spm('Pointer','Watch');

spm_slice_timing(P, Seq, refslice, timing, job.prefix);

spm('Pointer');

for i = 1:numel(job.scans)
    out(i).files = spm_file(job.scans{i}, 'prefix',job.prefix);
end

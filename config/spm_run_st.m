function out = spm_run_st(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_st.m 1185 2008-03-04 16:31:21Z volkmar $

job = varargin{1};
Seq = job.so;
TR  = job.tr;
TA  = job.ta;
nslices   = job.nslices;
refslice  = job.refslice;
timing(2) = TR - TA;
timing(1) = TA / (nslices -1);

for i = 1:length(job.scans)
    P   = strvcat(job.scans{i});
    spm_slice_timing(P,Seq,refslice,timing,job.prefix);
    for k = 1:numel(job.scans{i})
        [p n e v] = spm_fileparts(job.scans{i}{k});
        out(i).files{k} = fullfile(p, sprintf('%s%s%s%s', job.prefix, n, ...
                                              e, v));
    end;
end
return;

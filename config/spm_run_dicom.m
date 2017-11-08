function out = spm_run_dicom(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_dicom.m 7201 2017-11-08 11:13:25Z guillaume $


if ~isempty(job.outdir{1})
    out_dir = job.outdir{1};
else
    out_dir = pwd;
end

if job.convopts.icedims
    root_dir = ['ice' job.root];
else
    root_dir = job.root;
end

essentials = ~job.convopts.meta;
hdr = spm_dicom_headers(char(job.data),essentials);
sel = true(size(hdr));
if ~isempty(job.protfilter) && ~strcmp(job.protfilter, '.*')
    psel   = cellfun(@(h)isfield(h, 'ProtocolName'), hdr);
    ssel   = ~psel & cellfun(@(h)isfield(h, 'SequenceName'), hdr);
    pnames = cell(size(hdr));
    pnames(psel) = cellfun(@(h)subsref(h, substruct('.','ProtocolName')), hdr(psel), 'UniformOutput', false);
    pnames(ssel) = cellfun(@(h)subsref(h, substruct('.','SequenceName')), hdr(ssel), 'UniformOutput', false);
    sel(psel|ssel) = ~cellfun(@isempty,regexp(pnames(psel|ssel), job.protfilter));
end
out = spm_dicom_convert(hdr(sel),'all',root_dir,job.convopts.format,out_dir,job.convopts.meta);

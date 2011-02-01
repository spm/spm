function out = spm_run_preproc(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_preproc.m 4185 2011-02-01 18:46:18Z guillaume $

job.opts.tpm = char(job.opts.tpm);
if isfield(job.opts,'msk'),
    job.opts.msk = char(job.opts.msk);
end;
for i=1:numel(job.data),
    res           = spm_preproc(job.data{i},job.opts);
    [out(1).sn{i},out(1).isn{i}]   = spm_prep2sn(res);
    [pth,nam]     = spm_fileparts(job.data{i});
    out(1).snfile{i} = fullfile(pth,[nam '_seg_sn.mat']);
    savefields(out(1).snfile{i},out(1).sn{i});
    out(1).isnfile{i} = fullfile(pth,[nam '_seg_inv_sn.mat']);
    savefields(out(1).isnfile{i},out(1).isn{i});
end;
spm_preproc_write(cat(2,out.sn{:}),job.output);
% Guess filenames
opts  = job.output;
sopts = [opts.GM;opts.WM;opts.CSF];
for i=1:numel(job.data)
    [pth,nam,ext,num] = spm_fileparts(job.data{i});
    if opts.biascor,
        out(1).biascorr{i,1} = ...
            fullfile(pth, sprintf('m%s%s', nam, ext));
    end;
    for k1=1:3,
        if sopts(k1,3),
            out(1).(sprintf('c%d',k1)){i,1} = ...
                fullfile(pth, sprintf('c%d%s%s', k1, nam, ext));
        end;
        if sopts(k1,2),
            out(1).(sprintf('wc%d',k1)){i,1} = ...
                fullfile(pth, sprintf('wc%d%s%s', k1, nam, ext));
        end;
        if sopts(k1,1),
            out(1).(sprintf('mwc%d',k1)){i,1} = ...
                fullfile(pth, sprintf('mwc%d%s%s', k1, nam, ext));
        end;
    end;
end;
return;

%==========================================================================
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end
fn = fieldnames(p);
if numel(fn)==0, return; end
for i=1:length(fn)
    eval([fn{i} '= p.' fn{i} ';']);
end
if spm_check_version('matlab','7') >= 0
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end

return;

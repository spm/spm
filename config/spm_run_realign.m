function out = spm_run_realign(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_realign.m 7141 2017-07-26 09:05:05Z guillaume $


P = cell(size(job.data));
for i=1:numel(job.data)
    P{i} = char(job.data{i});
end

%-Realign
%--------------------------------------------------------------------------
if isfield(job,'eoptions')
    flags.quality = job.eoptions.quality;
    flags.fwhm    = job.eoptions.fwhm;
    flags.sep     = job.eoptions.sep;
    flags.rtm     = job.eoptions.rtm;
    flags.PW      = char(job.eoptions.weight);
    flags.interp  = job.eoptions.interp;
    flags.wrap    = job.eoptions.wrap;
    
    spm_realign(P, flags);
end

%-Reslice
%--------------------------------------------------------------------------
if isfield(job,'roptions')
    P            = char(P);
    flags.mask   = job.roptions.mask;
    flags.interp = job.roptions.interp;
    flags.which  = job.roptions.which;
    flags.wrap   = job.roptions.wrap;
    flags.prefix = job.roptions.prefix;
    
    spm_reslice(P, flags);
end

%-Dependencies
%--------------------------------------------------------------------------
if isempty(P), out = struct([]); return; end
if isfield(job,'eoptions')
    for i=1:numel(job.data)
        out.sess(i).cfiles = job.data{i};
        out.sess(i).rpfile{1} = spm_file(job.data{i}{1}, 'prefix','rp_', 'ext','.txt');
    end
end

if isfield(job,'roptions')
    if job.roptions.which(1) == 1, s = 1; else, s = 0; end
    if ischar(job.data{1}), job.data = {job.data}; end
    for k=1:numel(job.data)
        rfiles = cell(numel(job.data{k})-s,1);
        for i=1:numel(rfiles)
            rfiles{i} = spm_file(job.data{k}{i+s}, 'prefix',job.roptions.prefix);
        end
        if isfield(job,'eoptions')
            out.sess(k).rfiles = rfiles;
        else
            out.rfiles = rfiles;
        end
    end
    if job.roptions.which(2)
        out.rmean{1} = spm_file(job.data{1}{1}, 'prefix','mean', 'number','');
    end
end

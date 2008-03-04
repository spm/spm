function out = spm_run_realign_estimate(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_realign_estimate.m 1185 2008-03-04 16:31:21Z volkmar $

job           = varargin{1};
P             = {};
for i=1:length(job.data),
    P{i}  = strvcat(job.data{i});
end;
flags.quality = job.eoptions.quality;
flags.fwhm    = job.eoptions.fwhm;
flags.sep     = job.eoptions.sep;
flags.rtm     = job.eoptions.rtm;
flags.PW      = strvcat(job.eoptions.weight);
flags.interp  = job.eoptions.interp;
flags.wrap    = job.eoptions.wrap;
spm_realign(P,flags);
for i=1:numel(job.data)
    out.sess(i).cfiles = job.data{i};
    [pth,nam,ext,num]  = spm_fileparts(job.data{i}{1});
    out.sess(i).rpfile{1} =  fullfile(pth, sprintf('rp_%s.txt', nam));
end;
return;

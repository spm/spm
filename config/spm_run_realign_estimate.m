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

% $Id: spm_run_realign_estimate.m 4152 2011-01-11 14:13:35Z volkmar $

job           = varargin{1};
P             = cell(size(job.data));
for i=1:length(job.data),
    P{i}  = char(job.data{i});
end;
flags.quality = job.eoptions.quality;
flags.fwhm    = job.eoptions.fwhm;
flags.sep     = job.eoptions.sep;
flags.rtm     = job.eoptions.rtm;
flags.PW      = char(job.eoptions.weight);
flags.interp  = job.eoptions.interp;
flags.wrap    = job.eoptions.wrap;
spm_realign(P,flags);
for i=1:numel(job.data)
    out.sess(i).cfiles = job.data{i};
    [pth,nam,ext,num]  = spm_fileparts(job.data{i}{1});
    out.sess(i).rpfile{1} =  fullfile(pth, sprintf('rp_%s.txt', nam));
end;
return;

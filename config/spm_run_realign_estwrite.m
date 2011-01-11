function out = spm_run_realign_estwrite(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_realign_estwrite.m 4152 2011-01-11 14:13:35Z volkmar $

job           = varargin{1};
P             = cell(size(job.data));
for i=1:numel(job.data),
    P{i} = char(job.data{i});
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

P            = char(P);
flags.mask   = job.roptions.mask;
flags.interp = job.roptions.interp;
flags.which  = job.roptions.which;
flags.wrap   = job.roptions.wrap;
flags.prefix = job.roptions.prefix;
spm_reslice(P,flags);
switch job.roptions.which(1),
    case 1,
        for k=1:numel(job.data)
            out.sess(k).rfiles = cell(numel(job.data{k})-1,1);
            for i=1:length(out.sess(k).rfiles),
                [pth,nam,ext,num] = spm_fileparts(job.data{k}{i+1});
                out.sess(k).rfiles{i} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
            end;
        end;
    otherwise,
        for k=1:numel(job.data)
            out.sess(k).rfiles = cell(numel(job.data{k}),1);
            for i=1:length(out.sess(k).rfiles),
                [pth,nam,ext,num] = spm_fileparts(job.data{k}{i});
                out.sess(k).rfiles{i} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
            end;
        end;
end;
if job.roptions.which(2),
    [pth,nam,ext,num] = spm_fileparts(job.data{1}{1});
    out.rmean{1} = fullfile(pth,['mean', nam, ext, num]);
end;
return;

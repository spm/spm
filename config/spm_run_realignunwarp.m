function out = spm_run_realignunwarp(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Darren R. Gitelman
% $Id: spm_run_realignunwarp.m 4185 2011-02-01 18:46:18Z guillaume $

job = varargin{1};

% assemble flags
%-----------------------------------------------------------------------
% assemble realignment estimation flags.
flags.quality = job.eoptions.quality;
flags.fwhm    = job.eoptions.fwhm;
flags.sep     = job.eoptions.sep;
flags.rtm     = job.eoptions.rtm;
flags.PW      = char(job.eoptions.weight);
flags.interp  = job.eoptions.einterp;
flags.wrap    = job.eoptions.ewrap;

uweflags.order     = job.uweoptions.basfcn;
uweflags.regorder  = job.uweoptions.regorder;
uweflags.lambda    = job.uweoptions.lambda;
uweflags.jm        = job.uweoptions.jm;
uweflags.fot       = job.uweoptions.fot;

if ~isempty(job.uweoptions.sot)
    cnt = 1;
    for i=1:size(job.uweoptions.sot,2)
        for j=i:size(job.uweoptions.sot,2)
            sotmat(cnt,1) = job.uweoptions.sot(i);
            sotmat(cnt,2) = job.uweoptions.sot(j);
            cnt = cnt+1;
        end
    end
else
    sotmat = [];
end
uweflags.sot       = sotmat;
uweflags.fwhm      = job.uweoptions.uwfwhm;
uweflags.rem       = job.uweoptions.rem;
uweflags.noi       = job.uweoptions.noi;
uweflags.exp_round = job.uweoptions.expround;

uwrflags.interp    = job.uwroptions.rinterp;
uwrflags.wrap      = job.uwroptions.wrap;
uwrflags.mask      = job.uwroptions.mask;
uwrflags.which     = job.uwroptions.uwwhich(1);
uwrflags.mean      = job.uwroptions.uwwhich(2);
uwrflags.prefix    = job.uwroptions.prefix;

if uweflags.jm == 1
    uwrflags.udc = 2;
else
    uwrflags.udc = 1;
end
%-----------------------------------------------------------------------

% assemble files
%-----------------------------------------------------------------------
P   = cell(size(job.data));
sfP = cell(size(job.data));
for i = 1:numel(job.data)
    P{i} = char(job.data(i).scans{:});
    if ~isempty(job.data(i).pmscan)
        sfP{i} = job.data(i).pmscan{1};
    else
        sfP{i} = [];
    end
end

% realign
%-----------------------------------------------------------------------
spm_realign(P,flags);

for i = 1:numel(P)
    uweflags.sfP = sfP{i};

    % unwarp estimate
    %-------------------------------------------------------------------
    tmpP = spm_vol(P{i}(1,:));
    uweflags.M = tmpP.mat;
    ds = spm_uw_estimate(P{i},uweflags);
    out.sess(i).ds = ds;
    [path,name] = fileparts(P{i}(1,:));
    out.sess(i).dsfile{1} =  fullfile(path,[name '_uw.mat']);

    if spm_check_version('matlab','7') >= 0
        save(out.sess(i).dsfile{1},'-V6','ds');
    else
        save(out.sess(i).dsfile{1},'ds');
    end
end

% unwarp write - done at the single subject level since Batch
% forwards one subjects data at a time for analysis, assuming
% that subjects should be grouped as new spatial nodes. Sessions
% should be within subjects.
%-----------------------------------------------------------------------
spm_uw_apply(cat(2,out.sess.ds),uwrflags);
switch job.uwroptions.uwwhich(1)
    case 0
        out.sess.uwrfiles  = {};
    case 2
        for i = 1:numel(P)
            out.sess(i).uwrfiles = cell(size(P{i},1),1);
            for j=1:size(P{i},1)
                [pth,nam,ext,num] = spm_fileparts(deblank(P{i}(j,:)));
                out.sess(i).uwrfiles{j} = fullfile(pth,[job.uwroptions.prefix, ...
                    nam, ext, num]);
            end
        end
end
if job.uwroptions.uwwhich(2)
    [pth,nam,ext,num] = spm_fileparts(deblank(P{1}(1,:)));
    out.meanuwr{1} = fullfile(pth,['mean', job.uwroptions.prefix, nam, ext, num]);
end

return;

function out = spm_run_realignunwarp(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________

% Darren R. Gitelman
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


%-Assemble flags
%--------------------------------------------------------------------------
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

uwrflags.jm        = job.uwroptions.jm;
uwrflags.interp    = job.uwroptions.rinterp;
uwrflags.wrap      = job.uwroptions.wrap;
uwrflags.mask      = job.uwroptions.mask;
uwrflags.which     = job.uwroptions.uwwhich(1);
uwrflags.mean      = job.uwroptions.uwwhich(2);
uwrflags.prefix    = job.uwroptions.prefix;

%-Assemble files
%--------------------------------------------------------------------------
P   = cell(size(job.data));
sfP = cell(size(job.data));
for i = 1:numel(job.data)
    P{i} = char(job.data(i).scans{:});
    if ~isempty(job.data(i).pmscan)
        sfP{i} = job.data(i).pmscan{1};
    else
        sfP{i} = [];
        %fprintf('\n\n\n WARNING: Voxel displacement map (vdm*) not detected\n\n\n');
    end
end

%-Realign
%--------------------------------------------------------------------------
spm_realign(P,flags);

%-Unwarp Estimate
%--------------------------------------------------------------------------
M1  = spm_get_space(deblank(P{1}(1,:))); % Space of first image overall
for i = 1:numel(P)
    uweflags.sfP = sfP{i};
    P1           = deblank(P{i}(1,:));
    if isempty(spm_file(P1,'number'))
        P1 = spm_file(P1,'number',1);
    end
    VP1        = spm_vol(P1);
    uweflags.M = VP1.mat;

    % Move the static field map for this session to match the
    % first volume of the session
    if ~isempty(uweflags.sfP)
        uweflags.sfP     = spm_vol(uweflags.sfP);            % Map the file
        uweflags.sfP.mat = (uweflags.M/M1)*uweflags.sfP.mat; % Adjust its position
    end

    ds         = spm_uw_estimate(P{i},uweflags);
    sess(i).ds = ds;
    dsfile     = spm_file(P{i}(1,:), 'suffix','_uw', 'ext','.mat');
    save(dsfile,'ds', spm_get_defaults('mat.format'));
end

%-Unwarp Write - Sessions should be within subjects
%--------------------------------------------------------------------------
spm_uw_apply(cat(2,sess.ds),uwrflags);

%-Dependencies
%--------------------------------------------------------------------------
for i=1:numel(P)
    out.sess(i).rpfile{1} = spm_file(P{i}(1,:), 'prefix','rp_', 'ext','.txt');
    out.sess(i).dsfile{1} = spm_file(P{i}(1,:), 'suffix','_uw', 'ext','.mat');
end

switch job.uwroptions.uwwhich(1)
    case 0
        out.sess.uwrfiles  = {};
    case 2
        for i=1:numel(P)
            out.sess(i).uwrfiles = spm_file(cellstr(P{i}), ...
                'prefix',job.uwroptions.prefix);
        end
end
if job.uwroptions.uwwhich(2)
    out.meanuwr{1} = spm_file(P{1}(1,:), ...
        'prefix',['mean' job.uwroptions.prefix], ...
        'number','');
end

function out = cfg_run_runjobs(job)

% Initialise, fill, save and run a job with repeated inputs.
% To make use of possible parallel execution of independent jobs, all
% repeated jobs are filled first and (if successfully filled) run as one
% large job.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_runjobs.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

sts = true;
if isfield(job.save, 'savejobs')
    [p n e v] = fileparts(job.save.savejobs.outstub);
    outfmt = fullfile(job.save.savejobs.outdir{1}, sprintf('%s_%%0%dd.m', n, ceil(log10(numel(job.inputs))+1)));
end;
hjobs = {};
out.jobfiles = {};
for cr = 1:numel(job.inputs)
    cjob = cfg_util('initjob', job.jobs);
    for ci = 1:numel(job.inputs{cr})
        fn = fieldnames(job.inputs{cr}{ci});
        inp{ci} = job.inputs{cr}{ci}.(fn{1});
    end;
    sts1 = cfg_util('filljob', cjob, inp{:});
    if sts1
        [un hjobs{end+1}] = cfg_util('harvest', cjob);
    end;
    sts = sts && sts1;
    if isfield(job.save, 'savejobs')
        out.jobfiles{cr} = sprintf(outfmt, cr);
        cfg_util('savejob', cjob, out.jobfiles{cr});
    end;
    cfg_util('deljob', cjob);
end;
if sts || strcmp(job.missing,'skip')
    cjob = cfg_util('initjob', hjobs);
    cfg_util('run', cjob);
    if isfield(job.save, 'savejobs')
        [p n e v] = fileparts(job.save.savejobs.outstub);
        out.jobrun{1} = fullfile(p, [n '_run.m']);
        cfg_util('saverun', cjob, out.jobrun{1});
    end;
    cfg_util('deljob', cjob);
end;
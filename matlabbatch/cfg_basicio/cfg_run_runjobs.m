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
% $Id: cfg_run_runjobs.m 6514 2015-08-06 10:07:53Z volkmar $

rev = '$Rev: 6514 $'; %#ok

if isfield(job.save, 'savejobs')
    [p, n, e] = fileparts(job.save.savejobs.outstub);
    outfmt = strrep(fullfile(job.save.savejobs.outdir{1}, sprintf('%s_%%0%dd.m', n, ceil(log10(numel(job.inputs))+1))), '\', '\\');
end;
hjobs = cell(size(job.inputs));
sts = false(size(job.inputs));
out.jobfiles = {};
for cr = 1:numel(job.inputs)
    cjob = cfg_util('initjob', job.jobs);
    inp = cell(size(job.inputs{cr}));
    for ci = 1:numel(job.inputs{cr})
        fn = fieldnames(job.inputs{cr}{ci});
        inp{ci} = job.inputs{cr}{ci}.(fn{1});
    end;
    sts(cr) = cfg_util('filljob', cjob, inp{:});
    if sts(cr)
        [un, hjobs{cr}] = cfg_util('harvest', cjob);
    end;
    if isfield(job.save, 'savejobs')
        out.jobfiles{cr} = sprintf(outfmt, cr);
        cfg_util('savejob', cjob, out.jobfiles{cr});
    end;
    cfg_util('deljob', cjob);
end;
if all(sts)
    % keep all hjobs;
elseif any(sts) && strcmp(job.missing,'skip')
    hjobs = hjobs(sts);
else
    hjobs = {};
end
if ~isempty(hjobs)
    cjob = cfg_util('initjob', hjobs);
    cfg_util('run', cjob);
    out.jout = cfg_util('getalloutputs', cjob);
%    if isfield(job.save, 'savejobs')
%        [p n e] = fileparts(job.save.savejobs.outstub);
%        out.jobrun{1} = fullfile(p, [n '_run.m']);
%        cfg_util('saverun', cjob, out.jobrun{1});
%    end;
    cfg_util('deljob', cjob);
end;
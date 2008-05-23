function dep = cfg_vout_named_dir(job)

% Define virtual outputs for cfg_run_named_dir. Dir names can either be
% assigned to a cfg_dirs input or to a evaluated cfg_entry. Dir indices
% can be assigned to any numeric or evaluated cfg_entry item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_named_dir.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

if strcmp(job.name,'<UNDEFINED>') || isempty(job.name) || isa(job.name, 'cfg_dep')
    dirname = 'Directory';
else
    dirname = job.name;
end;

for k = 1:numel(job.dirs)
    dep(k) = cfg_dep;
    dep(k).sname = sprintf('%s(%d)', dirname, k);
    dep(k).src_output = substruct('.','dirs','{}',{k});
    dep(k).tgt_spec   = cfg_findspec({{'filter','dir','strtype','e'}});
end;

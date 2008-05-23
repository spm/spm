function dep = cfg_vout_save_vars(job)

% Define virtual output for cfg_run_save_vars. Output can be passed on to
% either a cfg_file or an evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_save_vars.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

if strcmp(job.name,'<UNDEFINED>') || isempty(job.name) || isa(job.name, 'cfg_dep')
    savename = '';
else
    savename = job.name;
end;

dep = cfg_dep;
dep.sname = sprintf('Save Variables: %s', savename);
dep.src_output = substruct('.','file');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
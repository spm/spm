function dep = cfg_vout_dir_move(job)

% Define virtual output for cfg_run_dir_move. Output can be passed on to
% either a cfg_files or an evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_dir_move.m 6920 2016-11-02 14:48:01Z guillaume $

rev = '$Rev: 6920 $'; %#ok

if ~isfield(job.action,'delete')
    dep = cfg_dep;
    dep.sname = 'Moved/Copied Directory';
    dep.src_output = substruct('.','dir');
    dep.tgt_spec   = cfg_findspec({{'filter','dir','strtype','e'}});
else
    dep = [];
end
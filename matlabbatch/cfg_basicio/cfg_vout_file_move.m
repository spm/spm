function dep = cfg_vout_move_file(job)

% Define virtual output for cfg_run_move_file. Output can be passed on to
% either a cfg_files or an evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_file_move.m 1456 2008-04-21 15:03:41Z volkmar $

rev = '$Rev: 1456 $';

if isfield(job.action,'moveto')
    dep = cfg_dep;
    dep.sname = 'Moved Files';
    dep.src_output = substruct('.','files');
    dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
else
    dep = [];
end;
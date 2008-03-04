function dep = cfg_vout_mkdir(job)

% Define virtual outputs for cfg_run_mkdir. The directory name can either
% be assigned to a cfg_files directory input or to a evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_mkdir.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';
dep            = cfg_dep;
if isa(job.name, 'cfg_dep')
    dep.sname  = sprintf('Make Directory');
else
    dep.sname  = sprintf('Make Directory ''%s''', job.name);
end;
dep.src_output = substruct('.','dir');
dep.tgt_spec   = cfg_findspec({{'filter','dir', 'strtype','e'}});

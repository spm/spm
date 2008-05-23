function dep = cfg_vout_file_fplist(job)

% function dep = cfg_vout_file_fplist(job)
%
% Virtual outputs for cfg_run_file_fplist. Struct with fields .files and
% .dirs. See help on cfg_getfile.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_file_fplist.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

dep(1)            = cfg_dep;
dep(1).sname      = 'Selected Files';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
dep(2)            = cfg_dep;
dep(2).sname      = 'Subdirectories';
dep(2).src_output = substruct('.','dirs');
dep(2).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

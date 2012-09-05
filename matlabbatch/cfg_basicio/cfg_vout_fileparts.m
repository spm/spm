function dep = cfg_vout_fileparts(job)

% Define virtual outputs for cfg_run_fileparts. The path names can either
% be assigned to a cfg_files input or to an evaluated cfg_entry. File names
% and extensions can only be assigned to an evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_fileparts.m 4899 2012-09-05 13:44:17Z volkmar $

rev = '$Rev: 4899 $'; %#ok
dep(1)            = cfg_dep;
dep(1).sname      = sprintf('Directories');
dep(1).src_output = substruct('.','p');
dep(1).tgt_spec   = cfg_findspec({{'class','cfg_files', 'strtype','e'}});
dep(2)            = cfg_dep;
dep(2).sname      = sprintf('Directories (unique)');
dep(2).src_output = substruct('.','up');
dep(2).tgt_spec   = cfg_findspec({{'class','cfg_files', 'strtype','e'}});
dep(3)            = cfg_dep;
dep(3).sname      = sprintf('Filenames');
dep(3).src_output = substruct('.','n');
dep(3).tgt_spec   = cfg_findspec({{'strtype','e'}});
dep(4)            = cfg_dep;
dep(4).sname      = sprintf('Extensions');
dep(4).src_output = substruct('.','e');
dep(4).tgt_spec   = cfg_findspec({{'strtype','e'}});

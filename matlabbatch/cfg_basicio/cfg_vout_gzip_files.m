function dep = cfg_vout_gzip_files(job)

% Define virtual outputs for "Gzip Files". File names can either be
% assigned to a cfg_files input or to a evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_gzip_files.m 6627 2015-12-04 08:20:49Z volkmar $

rev = '$Rev: 6627 $'; %#ok

dep            = cfg_dep;
dep.sname      = 'Gzipped Files';
dep.src_output = substruct('()',{':'});
dep.tgt_spec   = cfg_findspec({{'class','cfg_files', 'strtype','e'}});

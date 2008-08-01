function dep = cfg_vout_file_filter(job)

% Define virtual outputs for cfg_vout_file_filter. The file names can either
% be assigned to a cfg_files input or to a evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_file_filter.m 1973 2008-08-01 11:52:41Z volkmar $

rev = '$Rev: 1973 $'; %#ok
dep            = cfg_dep;
dep.sname      = sprintf('Filtered Files');
dep.src_output = substruct('.','files');
if ischar(job.typ) && any(strcmpi(job.typ, {'any','image','nifti','xml','mat','batch','dir'}))
    dep.tgt_spec   = cfg_findspec({{'filter',job.typ, 'strtype','e'}});
else
    dep.tgt_spec   = cfg_findspec({{'filter','any', 'strtype','e'}});
end
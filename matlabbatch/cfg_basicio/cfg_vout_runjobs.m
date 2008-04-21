function dep = cfg_vout_runjobs(job)

% Return dependency to jobfiles, if files are to be saved.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_runjobs.m 1456 2008-04-21 15:03:41Z volkmar $

rev = '$Rev: 1456 $';

dep = cfg_dep;
if isfield(job.save, 'savejobs')
    dep.sname = 'Job Files';
    dep.src_output = substruct('.','outfiles');
    dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
else
    dep = dep(false);
end;
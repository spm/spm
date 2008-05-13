function dep = cfg_vout_runjobs(job)

% Return dependency to jobfiles, if files are to be saved.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_runjobs.m 1606 2008-05-13 06:07:01Z volkmar $

rev = '$Rev: 1606 $';

dep = cfg_dep;
if isfield(job.save, 'savejobs')
    dep(1).sname = 'Job Files';
    dep(1).src_output = substruct('.','jobfiles');
    dep(1).tgt_spec   = cfg_findspec({{'filter','batch','strtype','e'}});
    dep(2)       = cfg_dep;
    dep(2).sname = 'Resulting Job after Computation';
    dep(2).src_output = substruct('.','jobrun');
    dep(2).tgt_spec   = cfg_findspec({{'filter','batch','strtype','e'}});
else
    dep = dep(false);
end;
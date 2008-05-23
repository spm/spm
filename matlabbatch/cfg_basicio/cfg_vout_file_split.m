function dep = cfg_vout_file_split(job)

% Define virtual outputs for cfg_run_file_split. File names can either be
% assigned to a cfg_files input or to a evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_file_split.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

if strcmp(job.name,'<UNDEFINED>') || isempty(job.name) || isa(job.name, 'cfg_dep')
    setname = 'File set';
else
    setname = job.name;
end;
    
for k = 1:numel(job.index)
    dep(k) = cfg_dep;
    dep(k).sname = sprintf('%s (%d)', setname, k);
    dep(k).src_output = substruct('{}', {k});
    dep(k).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
end;
dep(k+1) = cfg_dep;
dep(k+1).sname = sprintf('%s (rem)', setname);
dep(k+1).src_output = substruct('{}', {k+1});
dep(k+1).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

function dep = cfg_vout_named_file(job)

% Define virtual outputs for cfg_run_named_file. File names can either be
% assigned to a cfg_files input or to a evaluated cfg_entry. File indices
% can be assigned to any numeric or evaluated cfg_entry item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_vout_named_file.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

if strcmp(job.name,'<UNDEFINED>') || isempty(job.name) || isa(job.name, 'cfg_dep')
    setname = 'File set';
else
    setname = job.name;
end;

nf = numel(job.files);
for k = 1:nf
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('%s(%d) - Files', setname, k);
    dep(k).src_output = substruct('.','files','{}',{k});
    dep(k).tgt_spec   = cfg_findspec({{'class','cfg_files', 'strtype','e'}});
end;
for k = 1:nf
    dep(nf+k)            = cfg_dep;
    dep(nf+k).sname      = sprintf('%s(%d) - File index', setname, k);
    dep(nf+k).src_output = substruct('.','index','{}',{k});
    dep(nf+k).tgt_spec   = cfg_findspec({{'strtype','n', 'strtype','w', ...
                        'strtype','i', 'strtype','r', 'strtype','e'}});
end;
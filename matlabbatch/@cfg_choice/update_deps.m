function item = update_deps(item, varargin)

% function item = update_deps(item, varargin)
% This function will run cfg_dep/update_deps in all leaf (cfg_entry,
% cfg_menu, cfg_files) nodes of a configuration tree and update their
% dependency information (mod_job_ids) if necessary.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: update_deps.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

val = subsref(item, substruct('.','val'));
for k = 1:numel(val)
    val{k} = update_deps(val{k}, varargin{:});
end;
item = subsasgn(item, substruct('.','val'), val);
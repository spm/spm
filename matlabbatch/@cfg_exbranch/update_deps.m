function item = update_deps(item, oid, nid)

% function item = update_deps(item, varargin)
% This function will run cfg_dep/update_deps in all leaf (cfg_entry,
% cfg_menu, cfg_files) nodes of a configuration tree and update their
% dependency information (mod_job_ids) if necessary. It will also update
% cfg_exbranch mod_job_ids (item.id and ids in item.tdeps and item.sdeps).
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: update_deps.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

id = subsref(item, substruct('.', 'id'));
for k = 1:numel(oid)
    if isequal(id, oid{k})
        item = subsasgn(item, substruct('.', 'id'), nid{k});
        for cs = 1:numel(item.sout)
            item.sout(cs).src_exbranch = nid{k};
        end;
        break;
    end;
end;
if ~isempty(item.tdeps)
    item.tdeps = update_deps(item.tdeps, oid, nid);
end;
if ~isempty(item.sdeps)
    item.sdeps = update_deps(item.sdeps, oid, nid);
end;

val = subsref(item, substruct('.','val'));
for k = 1:numel(val)
    val{k} = update_deps(val{k}, oid, nid);
end;
item = subsasgn(item, substruct('.','val'), val);
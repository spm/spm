function item = initialise(item, val, dflag)

% function item = initialise(item, val, dflag)
% Initialise a configuration tree with values. If val is a job
% struct/cell, only the parts of the configuration that are present in
% this job will be initialised.
% If val has the special value '<DEFAULTS>', the entire configuration
% will be updated with values from .def fields. If a .def field is
% present in a cfg_leaf item, the current default value will be inserted,
% possibly replacing a previously entered (default) value.
% dflag is ignored in a cfg_branch.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: initialise.m 2101 2008-09-16 13:56:26Z volkmar $

rev = '$Rev: 2101 $'; %#ok

if strcmp(val,'<DEFAULTS>')
    item = initialise_def(item, val, dflag);
elseif isstruct(val)
    item = initialise_job(item, val, dflag);
else
    cfg_message('matlabbatch:initialise', ...
                'Can not initialise %s: job is not a struct.', ...
                gettag(item));
end;

function item = initialise_def(item, val, dflag)
citem = subsref(item, substruct('.','val'));
for k = 1:numel(citem)
    citem{k} = initialise(citem{k}, val, dflag);
end;
item = subsasgn(item, substruct('.','val'), citem);

function item = initialise_job(item, val, dflag)
% Determine possible tags
vtags = fieldnames(val);

for k = 1:numel(item.cfg_item.val)
    % find field in val that corresponds to one of the branch vals
    vi = strcmp(gettag(item.cfg_item.val{k}), vtags);
    if any(vi) % field names are unique, so there will be at most one match
        item.cfg_item.val{k} = initialise(item.cfg_item.val{k}, ...
            val.(vtags{vi}), dflag);
    end;
end;


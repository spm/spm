function cc = cfg_struct2cfg(co, indent)

% Import a config structure into a matlabbatch class tree. Input structures
% are those generated from the configuration editor, cfg2struct methods or
% spm_jobman config structures.
%
% The layout of the configuration tree and the types of configuration items
% have been kept compatible to a configuration system and job manager
% implementation in SPM5 (Statistical Parametric Mapping, Copyright (C)
% 2005 Wellcome Department of Imaging Neuroscience). This code has been
% completely rewritten based on an object oriented model of the
% configuration tree.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_struct2cfg.m 1606 2008-05-13 06:07:01Z volkmar $

rev = '$Rev: 1606 $';

if nargin < 2
    indent = '';
end;

% set warning backtrace to off, save old state
wsts = warning('query','backtrace');
warning('backtrace','off');
% set warning status according to cfg_get_defaults
warning(cfg_get_defaults(['warning.' mfilename '.verb']), ...
    'matlabbatch:cfg_struct2cfg:verb');

%% Class of node
% Usually, the class is determined by the node type. Only for branches
% there is a distinction necessary between executable and non-executable
% branches in spm_jobman config files.

if strcmp(co.type, 'branch') && isfield(co, 'prog')
    typ = 'cfg_exbranch';
else
    if numel(co.type > 4) && strcmp(co.type(1:4), 'cfg_')
        typ = co.type;
    else
        typ = sprintf('cfg_%s', co.type);
    end;
end;

try
    warning('matlabbatch:cfg_struct2cfg:verb', ...
            '%sNode %s (%s): %s > %s', indent, co.tag, co.name, co.type, typ);
catch
    warning('matlabbatch:cfg_struct2cfg:verb', ...
            '%sNode UNKNOWN: %s > %s', indent, co.type, typ);
end;
eval(sprintf('cc = %s;', typ));

%% Import children
% for branches, repeats, choices children are collected first, before the
% object is created.

val = {};
values = {};
switch typ
    case {'cfg_branch','cfg_exbranch'}
        for k = 1:numel(co.val)
            val{k} = cfg_struct2cfg(co.val{k}, [indent ' ']);
        end;
        co.val = val;
    case {'cfg_repeat', 'cfg_choice'}
        if isfield(co, 'val')
            for k = 1:numel(co.val)
                val{k} = cfg_struct2cfg(co.val{k}, [indent ' ']);
            end;            
            co.val = val;
        end;
        for k = 1:numel(co.values)
            values{k} = cfg_struct2cfg(co.values{k}, [indent ' ']);
        end;
        co.values = values;
end;

%% Assign fields
% try to assign fields, give warnings if something goes wrong

co = rmfield(co, 'type');
fn = fieldnames(co);
% omit id field, it has a different meaning in cfg_items than in spm_jobman
% and it does not contain necessary information.
idind = strcmp('id',fn);
fn = fn(~idind);

% treat name and tag fields first
try
    cc.name = co.name;
    fn = fn(~strcmp('name', fn));
end;
try
    cc.tag = co.tag;
    fn = fn(~strcmp('tag', fn));
end;
% if present, treat num field before value assignments
nind = strcmp('num',fn);
if any(nind)
    if numel(co.num) == 1
        onum = co.num;
        if isfinite(co.num) && co.num > 0
            co.num = [co.num co.num];
        else
            co.num = [0 Inf];
        end;
        warning('matlabbatch:cfg_struct2cfg:verb', ...
                '(WW) Node %s / ''%s'' field num [%d] padded to %s', ...
                cc.tag, cc.name, onum, mat2str(co.num));
    end;
    cc = try_assign(cc,co,'num');
    % remove num field from list
    fn = fn(~nind);
end;
    
% if present, convert .def field
nind = strcmp('def',fn);
if any(nind)
    if ~isempty(co.def) && ischar(co.def)
        % assume SPM5 style defaults
        co.def = {@spm_get_defaults co.def};
    end;
    cc = try_assign(cc,co,'def');
    % remove def field from list
    fn = fn(~nind);
end;

for k = 1:numel(fn)
    cc = try_assign(cc,co,fn{k});
end;

warning('backtrace',wsts.state);

function cc = try_assign(cc, co, fn)
try
    cc.(fn) = co.(fn);
catch
    warning('matlabbatch:cfg_struct2cfg:verb', ...
            '(WW) Node %s / ''%s'' field %s: import failed.', cc.tag, ...
            cc.name, fn);
    tmp = textscan(evalc('disp(lasterr)'), '%s', 'delimiter', '\n');
    for l = 1:numel(tmp{1})
        warning('matlabbatch:cfg_struct2cfg:verb', ...
                '(WW)       %s\n', tmp{1}{l});
    end;
end;

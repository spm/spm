function dep = subsasgn(dep, subs, varargin)

% function dep = subsasgn(dep, subs, varargin)
% subscript references we have to deal with are:
% one level
% dep.(field)   - i.e. struct('type',{'.'} ,'subs',{field})
% dep(idx)      - i.e. struct('type',{'()'},'subs',{idx})
% two levels
% dep(idx).(field)
%
% to be dealt with elsewhere
% dep.(field){fidx}
% three levels
% dep(idx).(field){fidx}
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn.m 4863 2012-08-27 08:09:23Z volkmar $

rev = '$Rev: 4863 $'; %#ok

persistent my_cfg_dep
sflag = false;
if strcmpi(subs(1).type, '()')
    % check array boundaries, extend if necessary
    if (numel(subs(1).subs) == 1 && numel(dep) < max(subs(1).subs{1})) || ...
            (numel(subs(1).subs) > ndims(dep)) || ...
            (numel(subs(1).subs) == ndims(dep) && any(cellfun(@max, subs(1).subs) > size(dep)))
        if isempty(my_cfg_dep)
            my_cfg_dep = cfg_dep;
        end
        if isempty(dep)
            dep = my_cfg_dep;
        end
        [dep(subs(1).subs{:})] = deal(my_cfg_dep);
    end
    if numel(subs) == 1
        [dep(subs(1).subs{:})] = deal(varargin{:});
        return;
    else
        % select referenced objects from input array, run subsasgn on them an
        % put results back into input array
        odep  = dep;
        osubs = subs(1);
        dep   = dep(subs(1).subs{:});
        subs  = subs(2:end);
        sflag = true;
    end
end
if  strcmpi(subs(1).type, '.')
    % field assignment
    if numel(subs) > 1
        % Only part of field value(s) assigned. Get old field value(s),
        % assign values to it.
        val = {dep.(subs(1).subs)};
        if nargin == 3
            % Same value to assign to all field values.
            val = cellfun(@(cval)subsasgn(cval,subs(2:end),varargin{1}), val, 'UniformOutput', false);
        else
            % Different values to assign to each field value.
            val = cellfun(@(cval,cin)subsasgn(cval,subs(2:end),cin), val, varargin, 'UniformOutput', false);
        end
    else
        % Field values to assign
        val = varargin;
    end
    % Assign to fields
    [ok val] = valcheck(subs(1).subs, val);
    if ok
        if isempty(dep)
            dep = repmat(cfg_dep,size(val));
        end
        [dep.(subs(1).subs)] = deal(val{:});
    end
else
    cfg_message('matlabbatch:subsref:unknowntype', 'Bad subscript type ''%s''.',subs(1).type);
end
if sflag
    odep(osubs.subs{:}) = dep;
    dep = odep;
end

function [ok val] = valcheck(subs, val)
persistent local_mysubs_fields;
if ~iscell(local_mysubs_fields)
    local_mysubs_fields = mysubs_fields;
end
ok = true;
switch subs
    case {'tname','sname'}
        if ~all(cellfun(@ischar,val))
            cfg_message('matlabbatch:subsasgn:name', 'Value for field ''%s'' must be a string.', ...
                subs);
            ok = false;
        end
    case {'tgt_spec'}
        sel = cellfun(@isempty,val);
        if any(sel)
            [val{sel}] = deal(cfg_findspec);
        end
        ok = all(cellfun(@is_cfg_findspec,val));
    case local_mysubs_fields,
        sel = cellfun(@isempty,val);
        if any(sel)
            [val{sel}] = deal(struct('type',{}, 'subs',{}));
        end
        ok = all(cellfun(@(cval)(isstruct(cval) && all(isfield(cval,{'type','subs'}))),val));
        if ~ok
            cfg_message('matlabbatch:subsasgn:subs', ['Value for field ''%s'' must be a struct with' ...
                ' fields ''type'' and ''subs''.'], subs);
        end
    otherwise
        cfg_message('matlabbatch:subsasgn:unknownfield', 'Reference to unknown field ''%s''.', subs);
end

function ok = is_cfg_findspec(fs)
ok = iscell(fs) && ...
    all(cellfun(@(cs)(isstruct(cs) && ...
    ((numel(fieldnames(cs))==1 && any(isfield(cs,{'name','value'}))) || ...
    (numel(fieldnames(cs))==2 && all(isfield(cs,{'name','value'}))))),fs));
if ~ok
    cfg_message('matlabbatch:ok_subsasgn:tgt_spec', ...
        'Target specification must be a cfg_findspec.');
end
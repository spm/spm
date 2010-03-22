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
% $Id: subsasgn.m 3792 2010-03-22 13:11:36Z volkmar $

rev = '$Rev: 3792 $'; %#ok

%% One-level subscripts
%--------------------------------------------------------------------------
if numel(subs) == 1
    switch subs(1).type,
        case {'.'},
            if numel(dep) == nargin-2
                for k = 1:numel(dep)
                    % input checks
                    ok = true;
                    switch subs(1).subs
                        case {'tname','sname'}
                            if ~ischar(varargin{k})
                                cfg_message('matlabbatch:subsasgn:name', 'Value for field ''%s'' must be a string.', ...
                                    subs(1).subs);
                                ok = false;
                            end
                        case {'tgt_spec'}
                            if isempty(varargin{k})
                                varargin{k} = cfg_findspec;
                            else
                                ok = iscell(varargin{k});
                                if ok
                                    for l = 1:numel(varargin{k})
                                        ok = isstruct(varargin{k}{l}) && numel(fieldnames(varargin{k}{l}))==2 ...
                                              && all(isfield(varargin{k}{l},{'name', ...
                                                            'value'}));
                                        if ~ok
                                            break;
                                        end
                                    end
                                end
                                if ~ok
                                    cfg_message('matlabbatch:ok_subsasgn:tgt_spec', ...
                                            'Target specification must be a cfg_findspec.');
                                end
                            end
                        case subs_fields(dep),
                            if isempty(varargin{k})
                                varargin{k} = struct('type',{}, 'subs',{});
                            elseif ~(isstruct(varargin{k}) && isfield(varargin{k},'type') && isfield(varargin{k},'subs'))
                                cfg_message('matlabbatch:subsasgn:subs', ['Value for field ''%s'' must be a struct with' ...
                                    ' fields ''type'' and ''subs''.'], subs(1).subs);
                                ok = false;
                            end
                        otherwise
                            cfg_message('matlabbatch:subsasgn:unknownfield', 'Reference to unknown field ''%s''.', subs(1).subs);
                    end
                    if ok,
                        dep(k).(subs(1).subs) = varargin{k};
                    end
                end
            else
                cfg_message('matlabbatch:subsasgn:numel', 'In an assignment  A.X = B, the number of elements in A and B must be the same.');
            end
        case {'()'},
            if isempty(dep)
                dep = cfg_dep;
            end
            dep(subs(1).subs{:}) = varargin{:};
        case {'{}'},
            cfg_message('matlabbatch:subsasgn:notcell', 'Cell content reference from non cell-array object.');
        otherwise
            cfg_message('matlabbatch:subsref:unknowntype', 'Unknown subsref type: ''%s''. This should not happen.', subs(1).type);
    end
    return;
end

%% Canonicalise and check multi-level subscripts
%--------------------------------------------------------------------------
if numel(subs) >= 2 && strcmp(subs(1).type,'.')
    % Canonicalise subs, so that there are two levels within this class
    subs = [struct('type','()','subs',{{':'}}) subs];
end

if ~strcmp(subs(1).type,'()')
    cfg_message('matlabbatch:subsasgn:wrongsubs', 'Wrong subscript reference.');
end

item1 = subsref(dep, subs(1));
if nargin-2 ~= 1 && nargin-2 ~= numel(item1)
    cfg_message('matlabbatch:subsasgn:numel', 'In an assignment  A.X = B, the number of elements in A and B must be the same.');
end

%% Two- and multi-level subscripts
%--------------------------------------------------------------------------
for k = 1:numel(item1)
    if numel(subs) == 2
        if nargin-2 == 1
            item1(k) = subsasgn(item1(k),subs(2),varargin{1});
        else
            item1(k) = subsasgn(item1(k),subs(2),varargin{k});
        end
    else
        val = subsref(item1(k), subs(2));
        % explicitly use builtin subsasgn for subscripts into field values
        if nargin-2 == 1
            val = builtin('subsasgn',val,subs(3:end),varargin{1});
        else
            val = builtin('subsasgn',val,subs(3:end),varargin{k});
        end
        item1(k) = subsasgn(item1(k), subs(2), val);
    end
end
dep = subsasgn(dep, subs(1), item1);
return;

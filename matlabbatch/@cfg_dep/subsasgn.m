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
% $Id: subsasgn.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

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
                                warning('matlabbatch:cfg_dep:subsasgn:string', 'Value for field ''%s'' must be a string.', ...
                                    subs(1).subs);
                                ok = false;
                            end;
                        case {'tgt_spec'}
                            % do not check tgt_spec
                        case subs_fields(dep),
                            if ~(isstruct(varargin{k}) && isfield(varargin{k},'type') && isfield(varargin{k},'subs'))
                                warning('matlabbatch:cfg_dep:subsasgn:substruct', ['Value for field ''%s'' must be a struct with' ...
                                    ' fields ''type'' and ''subs''.'], subs(1).subs);
                                ok = false;
                            end;
                        otherwise
                            error('matlabbatch:subsasgn:unknownfield', 'Reference to unknown field ''%s''.', subs(1).subs);
                    end;
                    if ok,
                        dep(k).(subs(1).subs) = varargin{k};
                    end;
                end;
            else
                error('matlabbatch:subsasgn:numel', 'In an assignment  A.X = B, the number of elements in A and B must be the same.');
            end;
        case {'()'},
            if isempty(dep)
                dep = cfg_dep;
            end;
            dep(subs(1).subs{:}) = varargin{:};
        case {'{}'},
            error('matlabbatch:subsasgn:notcell', 'Cell content reference from non cell-array object.');
        otherwise
            error('matlabbatch:subsref:unknowntype', 'Unknown subsref type: ''%s''. This should not happen.', subs(1).type);
    end;
    return;
end;

%% Canonicalise and check multi-level subscripts
%--------------------------------------------------------------------------
if numel(subs) >= 2 && strcmp(subs(1).type,'.')
    % Canonicalise subs, so that there are two levels within this class
    subs = [struct('type','()','subs',{{':'}}) subs];
end;

if ~strcmp(subs(1).type,'()')
    error('matlabbatch:subsasgn:wrongsubs', 'Wrong subscript reference.');
end;

item1 = subsref(dep, subs(1));
if nargin-2 ~= 1 && nargin-2 ~= numel(item1)
    error('matlabbatch:subsasgn:numel', 'In an assignment  A.X = B, the number of elements in A and B must be the same.');
end;

%% Two- and multi-level subscripts
%--------------------------------------------------------------------------
for k = 1:numel(item1)
    if numel(subs) == 2
        if nargin-2 == 1
            item1(k) = subsasgn(item1(k),subs(2),varargin{1});
        else
            item1(k) = subsasgn(item1(k),subs(2),varargin{k});
        end;
    else
        val = subsref(item1(k), subs(2));
        % explicitly use builtin subsasgn for subscripts into field values
        if nargin-2 == 1
            val = builtin('subsasgn',val,subs(3:end),varargin{1});
        else
            val = builtin('subsasgn',val,subs(3:end),varargin{k});
        end;
        item1(k) = subsasgn(item1(k), subs(2), val);
    end;
end;
dep = subsasgn(dep, subs(1), item1);
return;

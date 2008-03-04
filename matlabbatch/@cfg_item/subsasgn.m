function item = subsasgn(item, subs, varargin)

% function item = subsasgn(item, subs, varargin)
% This function implements subsasgn for all classes derived from cfg_item.
% It relies on the capability of each class constructor to re-classify a
% struct object after a new value has been assigned to its underlying
% struct (This capability has to be implemented in the derived class).
% The structure of a configuration tree does not permit any arrays of
% cfg_item objects. Therefore, the only subscript reference and
% assignment within an cfg_item is a dot assignment to fields of this
% cfg_item. 
% Subscript references we have to deal with are:
% one level
% item.(field)   - i.e. struct('type',{'.'} ,'subs',{field})
%
% to be dealt with elsewhere
% item.(field){fidx}
% 
% In a future version, '()' and '{}' subscripts may be supported to
% access val fields of a cfg_item tree as if they were part of a
% harvested job. For cfg_branch objects (where dot assignments are used
% for val fields in their job tree) it is mandatory to index the job as a
% struct array to access harvested fields.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

if numel(item) ~= 1
    error('matlabbatch:subsasgn:numel', ...
          'Arrays of cfg_item objects not supported.');
end;

%% One-level subscripts
%--------------------------------------------------------------------------
if numel(subs) == 1
    switch subs(1).type,
        case {'.'},
            if numel(item) == 1 && nargin == 3
                % structify item to access fields of derived class
                sitem = struct(item);
                citem = class(item);
                switch citem
                    case 'cfg_exbranch',
                        par_class = 'cfg_branch';
                        pf1 = subs_fields(sitem.cfg_branch);
                        pf2 = subs_fields(cfg_item);
                        par_fields = {pf1{:} pf2{:}};
                    case 'cfg_item',
                        par_class = '';
                        par_fields = {};
                    otherwise
                        par_class = 'cfg_item';
                        par_fields = subs_fields(sitem.cfg_item);
                end;
                switch subs(1).subs
                    case par_fields,
                        % call subsasgn_check of derived class to check
                        % for correctness
                        [ok val] = subsasgn_check(item,subs,varargin{1});
                        if ok,
                            sitem.(par_class) = subsasgn(sitem.(par_class), subs, val);
                        end
                    case subs_fields(item),
                        [ok val] = subsasgn_check(item,subs,varargin{1});
                        if ok,
                            sitem.(subs(1).subs) = val;
                        end;
                    otherwise
                        error('matlabbatch:subsasgn:unknownfield', ...
                              ['Reference to unknown field ''%s''.\n' ...
                               'To assign to a field in the job structure, use a reference like ' ...
                               '''(x).%s''.'], subs(1).subs, subs(1).subs);
                end;
                % re-classify (the class constructor must support
                % struct inputs and re-classification)
                item = feval(citem,sitem);
            else
                error('matlabbatch:subsasgn:numel', ...
                      ['Array assignments not supported for fields of cfg_item objects.\n' ...
                       'To assign to a field in the job structure, use a reference like ' ...
                       '''(x).%s''.'], subs(1).subs);

            end;
        case {'()','{}'},
        warning('matlabbatch:subsref:jobsubs', ...
                'Subscript type ''%s'' reserved for future use.', subs(1).type);
        otherwise
            error('matlabbatch:subsref:unknowntype', ...
                  'Unknown subsref type: ''%s''. This should not happen.', subs(1).type);
    end;
    return;
end;

%% Two- and multi-level subscripts
%--------------------------------------------------------------------------
item = subsasgn_rec(item, subs, varargin{1});
return;

function it = subsasgn_rec(it, subs, val)
if numel(subs) == 1
    % final subscript
    try
        it = builtin('subsasgn', it, subs, val);
    catch
        % there may be an object method that is missed by builtin
        it = subsasgn(it, subs, val);
    end;
else
    try
        it1 = subsref(it, subs(1));
    catch
        % it does not yet contain contents at this subscript
        % builtin subsasgn will be clever enough to create the necessary
        % struct/fields that are referenced
        it = builtin('subsasgn', it, subs, val);
        return;
    end;
    it1 = subsasgn_rec(it1, subs(2:end), val);
    try
        % do we need this?
        it = builtin('subsasgn', it, subs(1), it1);
    catch
        % try to use class subsasgn
        it = subsasgn(it,subs(1), it1);
    end;
end;
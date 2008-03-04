function varargout = subsref(item, subs)

% function varargout = subsref(item, subs)
% subscript references we have to deal with are:
% one level
% item.(field)   - i.e. struct('type',{'.'} ,'subs',{field})
% item(idx)      - i.e. struct('type',{'()'},'subs',{idx})
% two levels
% item(idx).(field)
%
% to be dealt with elsewhere
% item.(field){fidx}
% three levels
% item(idx).(field){fidx}
% This function is identical for all classes derived from cfg_item, but it
% needs to be present in the class folder to access fields added by the
% derived class.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsref.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

switch subs(1).type,
    case {'.'},
        if numel(subs) > 1 && numel(item) > 1
            error('matlabbatch:subsref:multiref', 'Field reference for multiple structure elements that is followed by more reference blocks is an error.');
        end;
        switch subs(1).subs
            case subs_fields(item),
                for k = 1:numel(item)
                    val{k} = item(k).(subs(1).subs);
                end;
            case cat(2, subs_fields(item.cfg_branch), subs_fields(cfg_item)),
                for k = 1:numel(item)
                    val{k} = item(k).cfg_branch.(subs(1).subs);
                end;
            otherwise
                error('matlabbatch:subsref:unknownfield', ...
                      ['Reference to unknown field ''%s''.\nTo reference ' ...
                       'a field in the job structure, use a reference like ' ...
                       '''(x).%s'''], subs(1).subs, subs(1).subs);
        end;
    case {'()','{}'},
        val = subsref_job(item, subs, false);
    otherwise
        error('matlabbatch:subsref:unknowntype', ...
              'Unknown subsref type: ''%s''. This should not happen.', subs(1).type);
end
if strcmp(subs(1).type, '.') && numel(subs) > 1 
    % in this case, val has only one element, and subs(2:end) are indices into val{1}
    %    val{1} = builtin('subsref', val{1}, subs(2:end));
    % The line above does not seem to work, as MATLAB is not able to figure out
    % which subsref to call and when.
    for k = 2:numel(subs)
        val{1} = subsref(val{1}, subs(k));
    end;
end;

varargout = val;
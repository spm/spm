function varargout = subsref(dep, subs)

% function varargout = subsref(dep, subs)
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
% $Id: subsref.m 2512 2008-12-01 13:21:29Z volkmar $

rev = '$Rev: 2512 $'; %#ok

switch subs(1).type,
    case {'.'},
        if numel(subs) > 1 && numel(dep) > 1
            cfg_message('matlabbatch:subsref:multiref', 'Field reference for multiple structure elements that is followed by more reference blocks is an error.');
        end;
        switch subs(1).subs
            case subs_fields(dep),
                val = cell(size(dep));
                for k = 1:numel(dep)
                    val{k} = dep(k).(subs(1).subs);
                end;
            otherwise
                cfg_message('matlabbatch:subsref:unknownfield', 'Reference to unknown field ''%s''.', subs(1).subs);
        end;
    case {'()'},
        if numel(subs(1).subs) == 1 % vectorise output
            szi = numel(dep);
        else % sub-index output
            szi = size(dep);
        end;
        for k = 1:numel(szi)
            if ischar(subs(1).subs{k}) && strcmp(subs(1).subs{k},':')
                subs(1).subs{k} = 1:szi(k);
            end;
        end;        
        val = {dep(subs(1).subs{:})};
    case {'{}'}
        cfg_message('matlabbatch:subsref:notcell', 'Cell content reference from non cell-array object.');
    otherwise
        cfg_message('matlabbatch:subsref:unknowntype', 'Unknown subsref type: ''%s''. This should not happen.', subs(1).type);
end
if numel(subs) > 1 % in this case, val has only one element, and subs(2:end) are indices into val{1}
    val = {builtin('subsref', val{1}, subs(2:end))};
end;

varargout = val;
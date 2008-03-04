function sts = subsasgn_check_valcfg(subs,val,num)

% function sts = subsasgn_check_valcfg(subs,val,num)
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check_valcfg.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

sts = true;

if numel(subs) == 1
    % pass a cell of items
    if ~iscell(val)
        warning('matlabbatch:subsasgn_check_valcfg:iscell', 'Value must be a cell of ''cfg_item'' objects.');
        sts = false;
        return;
    end;
    % check whether number of items fits into num restrictions
    cn = numel(val);
    if (cn < num(1))||(cn > num(2))
        warning('matlabbatch:subsasgn_check_valcfg:numel', 'Number of values %d must be in range [%d %d].', cn, num(1), ...
            num(2));
        sts = false;
        return;
    end;
    % check whether each val is a menu config item
    for k=1:numel(val)
        if ~isa(val{k},'cfg_item')
            warning('matlabbatch:subsasgn_check_valcfg:cfg_item', 'Value for ''val{%d}'' must be a ''cfg_item''.',k);
            sts = false;
            return;
        end;
    end;

elseif numel(subs) == 2
    % pass either a cell of items or a single item
    if max(subs(2).subs{1}) > num(2)
        warning('matlabbatch:subsasgn_check_valcfg:numel', 'Number of values must be in range [%d %d].', num(1), ...
            num(2));
        sts = false;
        return;
    end;
    switch(subs(2).type)
        case {'()'}
            for k=1:numel(val)
                if ~isa(val{k},'cfg_item')
                    warning('matlabbatch:subsasgn_check_valcfg:cfg_item', 'Value for ''val{%d}'' must be a ''cfg_item''.', k);
                    sts = false;
                    return;
                end;
            end;
        case {'{}'}
            if ~isa(val,'cfg_item')
                warning('matlabbatch:subsasgn_check_valcfg:cfg_item', 'Value for ''val'' must be a ''cfg_item''.');
                sts = false;
                return;
            end;
    end;
end;

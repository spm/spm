function sts = subsasgn_check_valcfg(subs,val,num)

% function sts = subsasgn_check_valcfg(subs,val,num)
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check_valcfg.m 8190 2021-11-16 18:26:42Z guillaume $

rev = '$Rev: 8190 $'; %#ok

sts = true;

% pass a cell of items
if isa(val,'function_handle')
    val = feval(val);
end
if ~iscell(val)
    cfg_message('matlabbatch:checkval', ...
            'Value must be a cell of ''cfg_item'' objects.');
    sts = false;
    return;
end
% check whether number of items fits into num restrictions
cn = numel(val);
if (cn < num(1))||(cn > num(2))
    cfg_message('matlabbatch:checkval', ...
            'Number of values %d must be in range [%d %d].', cn, num(1), ...
        num(2));
    sts = false;
    return;
end
% check whether each val is a menu config item
for k=1:numel(val)
    if ~isa(val{k},'cfg_item')
        cfg_message('matlabbatch:checkval', ...
                'Value for ''val{%d}'' must be a ''cfg_item''.',k);
        sts = false;
        return;
    end
end
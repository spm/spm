function item = subsasgn_job(item, subs, val)

% function item = subsasgn_job(item, subs, val)
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. If subs is empty, then the subtree
% beginning at item will be initialised with val. Otherwise, subs(1)
% should have a subscript type of '.' in combination with a tagname from
% item.val. 
%
% This function is identical for cfg_branch and cfg_(m)choice classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_job.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

if isempty(subs)
    item = initialise(item, '<DEFAULTS>', false);
    item = initialise(item, val, false);
else
    [sts, vind] = checksubs_job(item, subs, false);
    if sts
        citem = subsasgn_job(subsref(item, substruct('.', treepart(item, false),'{}',{vind})), subs(2:end), val);
        item  = subsasgn(item, substruct('.', treepart(item, false), '{}', {vind}), citem);
    else
        [sts, vind] = checksubs_job(item, subs, true);
        if sts % should only be reached for cfg_(m)choice items
            citem = subsasgn_job(subsref(item, substruct('.', treepart(item, true),'{}',{vind})), subs(2:end), val);
            item  = subsasgn(item, substruct('.', treepart(item, false), '{}', {numel(item.cfg_item.val)+1}), citem);
        else
            cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
        end
    end
end

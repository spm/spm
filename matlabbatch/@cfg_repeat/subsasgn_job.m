function item = subsasgn_job(item, subs, val)

% function varargout = subsasgn_job(item, subs)
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. If subs is empty, then the subtree
% beginning at item will be initialised with val. Otherwise, subs(1)
% should have a subscript type of '()' or '{}', depending on the
% structure of the harvested job (cell or struct array). subs(2) should
% have a subscript type of '.' and the subscript reference should be a
% tagname from item.val.
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
elseif numel(subs) == 1
    if any(strcmp(subs(1).type, {'()', '{}'}))
        if max(subs(1).subs{1}) > numel(item.cfg_item.val)
            defsubs = setdiff((numel(item.cfg_item.val)+1):max(subs(1).subs{1}), subs(1).subs{1});
        else
            defsubs = [];
        end
        if strcmp(subs(1).type, '()')
            if numel(val) == 1
                val = repmat(val, 1, numel(subs(1).subs{1}));
            elseif numel(val) ~= numel(subs(1).subs{1})
                cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
            end
            if numel(item.values) == 1 && ~item.forcestruct
                % struct array reference or () subscript into cell array
                citem = initialise(item.values{1}, '<DEFAULTS>', false);
                item  = subsasgn(item, substruct('.', treepart(item, false), '()', {defsubs}), {citem});
                for k = 1:numel(subs(1).subs{1})
                    citem = initialise(item.values{1}, '<DEFAULTS>', false);
                    if isstruct(val)
                        val1 = val(k);
                    else
                        val1 = val{k};
                    end
                    citem = initialise(citem, val1, false);
                    item  = subsasgn(item, substruct('.', treepart(item, false), '{}', {subs(1).subs{1}(k)}), citem);                    
                end
            else
                % () subscript into cell array, no default items allowed
                if isempty(defsubs)
                    item1 = initialise(item, '<DEFAULTS>', false);
                    item1 = initialise(item1, val, false);
                    for k = 1:numel(subs(1).subs{1})
                        item = subsasgn(item, substruct('.', treepart(item, false), '{}', {subs(1).subs{1}(k)}), item1.cfg_item.val{k});
                    end
                else
                    cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
                end
            end
        elseif strcmp(subs(1).type, '{}')
            if numel(item.values) == 1 && ~item.forcestruct 
                if isa(item.values{1}, 'cfg_branch')
                    % struct array - no {} reference allowed
                    cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
                else
                    citem = initialise(item.values{1}, '<DEFAULTS>', false);
                    item  = subsasgn(item, substruct('.', treepart(item, false), '()', {defsubs}), {citem});
                    for k = 1:numel(subs(1).subs{1})
                        citem = initialise(item.values{1}, '<DEFAULTS>', false);
                        citem = initialise(citem, val, false);
                        item  = subsasgn(item, substruct('.', treepart(item, false), '{}', {subs(1).subs{1}(k)}), citem);
                    end                   
                end
            elseif isstruct(val)
                val1 = cell(size(val));
                for k = 1:numel(val)
                    val1{k} = val(k);
                end
                item1 = initialise(item, '<DEFAULTS>', false);
                item1 = initialise(item1, val1, false);
                for k = 1:numel(subs(1).subs{1})
                    item = subsasgn(item, substruct('.', treepart(item, false), '{}', {subs(1).subs{1}(k)}), item1.cfg_item.val{k});
                end
            else
                cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
            end
        end
    elseif strcmp(subs(1).type, '.') && (numel(item.values) == 1) && ~item.forcestruct && isa(item.values{1}, 'cfg_branch')
        if numel(item.cfg_item.val) == 0
            item.cfg_item.val{1} = initialise(item.values{1}, '<DEFAULTS>', false);
            item.cfg_item.val{1} = subsasgn_job(item.cfg_item.val{1}, subs, val);
        elseif numel(item.cfg_item.val) == 1
            item.cfg_item.val{1} = subsasgn_job(item.cfg_item.val{1}, subs, val);
        else
            cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
        end
    else
        cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
    end
else
    [sts, vind] = checksubs_job(item, subs, false);
    if (numel(item.values) == 1) && ~item.forcestruct && isa(item.values{1}, 'cfg_branch')
        nsubs = subs(2:end);
    else
        nsubs = subs(3:end);
    end
    if sts
        citem = subsasgn_job(subsref(item, substruct('.', treepart(item, false),'{}',{vind})), nsubs, val);
        item  = subsasgn(item, substruct('.', treepart(item, false), '{}', subs(1).subs), citem);
    else
        [sts, vind] = checksubs_job(item, subs, true);
        if sts
            if numel(item.values) == 1
                if max(subs(1).subs{1}) > numel(item.cfg_item.val)
                    defsubs = setdiff((numel(item.cfg_item.val)+1):max(subs(1).subs{1}), subs(1).subs{1});
                else
                    defsubs = [];
                end
                if numel(val) == 1
                    val = repmat(val, 1, numel(subs(1).subs{1}));
                elseif numel(val) ~= numel(subs(1).subs{1})
                    cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
                end
                citem = initialise(item.values{1}, '<DEFAULTS>', false);
                item  = subsasgn(item, substruct('.', treepart(item, false), '()', {defsubs}), {citem});
                citem = subsasgn_job(subsref(item, substruct('.', treepart(item, true),'{}',{vind})), nsubs, val);
                item  = subsasgn(item, substruct('.', treepart(item, false), '{}', subs(1).subs), citem);
            else
                if max(subs(1).subs{1}) > numel(item.cfg_item.val) + 1
                    cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
                else
                    citem = subsasgn_job(subsref(item, substruct('.', treepart(item, true),'{}',{vind})), nsubs, val);
                    item  = subsasgn(item, substruct('.', treepart(item, false), '{}', subs(1).subs), citem);
                end
            end
        else
            cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
        end
    end
end
function varargout = subsref_job(item, subs, c0)

% function [ritem varargout] = subsref_job(item, subs, c0)
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. If subs is empty, then the harvested subtree
% beginning at item will be returned. Otherwise, subs(1) should have a
% subscript type of '()' or '{}', depending on the structure of the
% harvested job (cell or struct array). If subs has more than one
% components, subs(2) should have a subscript type of '.' and the
% subscript reference should be a tagname from item.val.
% The third argument c0 is a copy of the entire job configuration. This
% is only used to reference dependencies properly.
% The first values returned is the referenced cfg_item object. The
% following values are the results of sub-referencing into item.val{1}.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsref_job.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

if isempty(subs)
    varargout{1} = item;
    [un, varargout{2}] = harvest(item, c0, false, false);
else
    [sts, vind] = checksubs_job(item, subs, false);
    if sts
        if numel(subs) == 1
            varargout{1} = item;
            [un, val] = harvest(item, c0, false, false);
            [varargout{2:nargout}] = cfg_callbuiltin('subsref', val, subs);
        else
            csubs = [substruct('.', treepart(item, false)) subs(1)];
            [varargout{1:nargout}] = subsref_job(subsref(item, csubs), subs(3:end), c0);
        end
    else
        [sts, vind] = checksubs_job(item, subs, true);
        if sts
            if numel(subs) == 1
                tn = tagnames(item, true);
                varargout{1} = item;
                varargout{2} = cell2struct(cell(size(tn)), tn, 2);
            else
                csubs = substruct('.', treepart(item, true), '{}', {vind});
                [varargout{1:nargout}] = subsref_job(subsref(item, csubs), subs(3:end), c0);
            end
        else
            cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));
        end
    end
end
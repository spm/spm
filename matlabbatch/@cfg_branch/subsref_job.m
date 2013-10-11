function varargout = subsref_job(item, subs, c0)

% function [ritem varargout] = subsref_job(item, subs, c0)
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. If subs is empty, then the harvested subtree
% beginning at item will be returned. Otherwise, subs(1) should have a
% subscript type of '.' in combination with a tagname from item.val.
% The third argument c0 is a copy of the entire job configuration. This
% is only used to reference dependencies properly.
% The first values returned is the referenced cfg_item object. The
% following values are the results of sub-referencing into item.val{x}.
%
% This function is identical for cfg_branch and cfg_(m)choice classes.
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
        csubs = substruct('.', treepart(item, false));
        citems = subsref(item, csubs);
        [varargout{1:nargout}] = subsref_job(citems{vind}, subs(2:end), c0);
    else
        [sts, vind] = checksubs_job(item, subs, true);
        if sts
            csubs = substruct('.', treepart(item, true));
            citems = subsref(item, csubs);
            [varargout{1:nargout}] = subsref_job(citems{vind}, subs(2:end), c0);
        else
            cfg_message('matlabbatch:subsref_job', 'Item ''%s'': invalid subscript reference.', gettag(item));        
        end
    end
end
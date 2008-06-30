function val = subsref_job(item, subs, dflag)

% function val = subsref_job(item, subs, dflag)
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. This function is only defined for in-tree nodes.
% The structure of a subscript reference is always a '()' or '{}'
% reference, followed by a '.' reference with the name of a child tag of
% the current node.
% If the referenced child is a leaf node, then its val{1} is
% returned. Else, if there are more subscripts, then subsref_job
% continues on the child node. Finally, if there are no more subscripts
% left, then the child cfg_item is returned.
% For cfg_choices, the sign of the first subscript index '()' or '{}' indicates
% whether the filled (ind(1) >= 0) or defaults (ind(1) < 0) part of the
% tree should be visited. Once a subscript points into the defaults tree,
% all subsequent subscript signs will be ignored and always the defaults
% tree will be visited.
% For cfg_branches, the first subscript '()' or '{}' is ignored.
% The second subscript '.' is evaluated against the tag names of the
% cfg_items in the .val or .values field, resp..
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsref_job.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

dflag = dflag || subs(1).subs{1} < 0;
tn = tagnames(item, dflag);
tp = treepart(item, dflag);
citems = subsref(item, substruct('.', tp));
if strcmp(subs(2).type, '.')
    cval = find(strcmp(subs(2).subs, tn));
    if ~isempty(cval)
        if isa(citems{cval},'cfg_intree')
            if numel(subs) > 2
                val = subsref_job(citems{cval}, subs(3:end), dflag);
            else
                val = citems{cval};
            end;
        else
            [tag val] = harvest(citems{cval},citems{cval}, false, false);
            if numel(subs) > 2 && ~strcmp(val1, '<UNDEFINED>')
                try
                    val = subsref(val1, subs(3:end));
                catch
                    cfg_message('matlabbatch:subsref', ...
                            'Subscript into value failed.');
                end;
            end;
        end;
    else
        cfg_message('matlabbatch:subsref', ...
              'Reference to unknown field ''%s'' in job structure.', ...
              subs(2).subs);
    end;
else
    cfg_message('matlabbatch:subsref', ...
          'Wrong subscript reference ''%s''.', subs(2).type);
end;
    
val = {val};
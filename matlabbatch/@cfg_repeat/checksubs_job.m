function [sts, vind] = checksubs_job(item, subs, dflag)

% function [sts vind] = checksubs_job(item, subs, dflag)
% Check whether a subscript reference is a valid reference in a job
% structure starting at item. subs(1) should have a subscript type of
% '()' or '{}', depending on the structure of the harvested job (cell or
% struct array). If subs has more than one components, subs(2) should
% have a subscript type of '.' and the subscript reference should be a
% tagname from item.val or item.values, depending on dflag.
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: checksubs_job.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

vind = [];
try
    sts = (((numel(item.values) == 1) && ~item.forcestruct && isa(item.values{1}, 'cfg_branch') && strcmp(subs(1).type, '()')) || any(strcmp(subs(1).type, {'()', '{}'}))) && (dflag || max(subs(1).subs{1}) <= numel(item.cfg_item.val));
    if sts && numel(subs) > 1
        if strcmp(subs(2).type, '.')
            if (numel(item.values) == 1) && ~item.forcestruct && isa(item.values{1}, 'cfg_branch') && strcmp(subs(1).type, '()')
                if dflag
                    vind = 1;
                else
                    vind = subs(1).subs{1};
                end
            else
                if dflag
                    vind = find(strcmp(subs(2).subs, tagnames(item, dflag)),1);
                    sts  = ~isempty(vind);
                else
                    sts = strcmp(subsref(tagnames(item, dflag), subs(1)), subs(2).subs);
                    if sts
                        vind = subs(1).subs{1};
                    end
                end
            end
        else
            sts = false;
            vind = [];
        end
    end
catch
    sts = false;
end
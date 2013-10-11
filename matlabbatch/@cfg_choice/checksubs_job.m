function [sts, vind] = checksubs_job(item, subs, dflag)

% function [sts vind] = checksubs_job(item, subs, dflag)
% Check whether a subscript reference is a valid reference in a job
% structure starting at item. subs(1) should have a subscript type of
% '.', and the subscript reference should be a tagname from item.val or
% item.values, depending on dflag.
%
% This function is identical for cfg_branch and cfg_(m)choice classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: checksubs_job.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

if strcmp(subs(1).type, '.')
    vind = find(strcmp(subs(1).subs, tagnames(item, dflag)),1);
    sts  = ~isempty(vind);
else
    sts = false;
    vind = [];
end

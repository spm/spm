function varargout = subsref_job(item, subs, c0)

% function [ritem varargout] = subsref_job(item, subs, c0)
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. This generic cfg_item method treats subs as a
% subscript reference into item.val{1}. This is suitable for all cfg_leaf
% items.
% The third argument c0 is a copy of the entire job configuration. This
% is only used to reference dependencies properly.
% The first value returned is the referenced cfg_item object. The
% following values are the results of sub-referencing into item.val{1}.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsref_job.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

varargout{1} = item;
[un val] = harvest(item, c0, false, false);
if isempty(subs)
    varargout{2} = val;
else
    [varargout{2:nargout}] = cfg_callbuiltin('subsref', val, subs);
end

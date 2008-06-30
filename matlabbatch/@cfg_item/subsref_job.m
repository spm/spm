function val = subsref_job(item, subs, dflag)

% function val = subsref_job(item, subs, dflag)
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. This function is only defined for in-tree nodes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsref_job.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

cfg_message('matlabbatch:subsref', ...
        'Subscript type ''%s'' reserved for future use.', subs(1).type);
val = {};
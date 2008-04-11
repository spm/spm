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
% $Id: subsref_job.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

warning('matlabbatch:subsref:jobsubs', ...
        'Subscript type ''%s'' reserved for future use.', subs(1).type);
val = {};
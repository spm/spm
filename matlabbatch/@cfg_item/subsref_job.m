function val = subsref_job(item, subs, dflag)

% function val = subsref_job(item, subs, dflag)
%
% Treat a subscript reference as a reference in a job structure instead
% of a cfg_item structure. This function is only defined for in-tree nodes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsref_job.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

warning('matlabbatch:subsref:jobsubs', ...
        'Subscript type ''%s'' reserved for future use.', subs(1).type);
val = {};
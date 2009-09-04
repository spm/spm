function cfg_run_assignin(job)

% Assign the value of job.output to a workspace variable job.name.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_assignin.m 3355 2009-09-04 09:37:35Z volkmar $

rev = '$Rev: 3355 $'; %#ok

% check for existence of variable
vars = evalin('base','feval(@who);');
% generate new name
name = genvarname(job.name, vars);
assignin('base', name, job.output);

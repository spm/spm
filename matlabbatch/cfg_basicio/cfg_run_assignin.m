function cfg_run_assignin(job)

% Assign the value of job.output to a workspace variable job.name.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_assignin.m 3434 2009-09-30 13:01:28Z volkmar $

rev = '$Rev: 3434 $'; %#ok

% check for existence of variable
vars = evalin('base','feval(@who);');
% generate new name
name = genvarname(job.name, vars);
if ~isequal(name,job.name)
    cfg_message('cfg_basicio:cfg_run_assignin:newname', ['Using ''%s'' ' ...
                        'instead of suggested variable name ''%s''.'], ...
                name, job.name);
end
assignin('base', name, job.output);

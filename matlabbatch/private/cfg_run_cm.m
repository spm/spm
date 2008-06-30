function cm = cfg_run_cm(cm, job)

% function cm = cfg_run_cm(cm, job)
% Run a module and return its output. Should really become a method of
% cfg_exbranch classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_cm.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

if isempty(cm.vout) && ~isempty(cm.vfiles);
    cfg_message('matlabbatch:deprecated:vfiles', ...
            'Using deprecated ''vfiles'' output in node ''%s''.', cm.tag);
    feval(cm.prog, job);
    cm.jout = struct('vfiles', {feval(cm.vfiles, ...
                                      job)});
elseif isempty(cm.sout)
    % no outputs specified
    feval(cm.prog, job);
    cm.jout = [];
else
    cm.jout = feval(cm.prog, job);
end;

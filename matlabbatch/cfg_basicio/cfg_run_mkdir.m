function out = cfg_run_mkdir(job)

% Make a directory and return its path in out.dir{1}.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_mkdir.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

out.dir{1} = fullfile(job.parent{1}, job.name);
if ~exist(out.dir{1}, 'dir')
    mkdir(out.dir{1});
end;
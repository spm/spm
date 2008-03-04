function cfg_run_cd(job)

% Make a directory and return its path in out.dir{1}.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_cd.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

cd(job.dir{1});
fprintf('Change Directory: New working directory\n\n     %s\n\n', ...
        job.dir{1})

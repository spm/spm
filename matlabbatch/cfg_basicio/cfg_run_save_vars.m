function out = cfg_run_save_vars(job)

% Save input variables to .mat file - either as a struct array, or as
% individual variables.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_save_vars.m 6929 2016-11-14 13:07:31Z guillaume $

rev = '$Rev: 6929 $'; %#ok

[p, n, e] = fileparts(job.name);
out.file{1} = fullfile(job.outdir{1}, [n '.mat']);
for k = 1:numel(job.vars)
    svar.(job.vars(k).vname) = job.vars(k).vcont; %#ok<STRNU>
end;
if job.saveasstruct
    save(out.file{1}, 'svar','-v6');
else
    save(out.file{1}, '-struct','svar','-v6');
end;

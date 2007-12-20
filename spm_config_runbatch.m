function opts = spm_config_runbatch
% Configuration file for running batched jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Darren Gitelman
% $Id: spm_config_runbatch.m 1032 2007-12-20 14:45:55Z john $

data.type = 'files';
data.name = 'Batch Files';
data.tag  = 'jobs';
data.filter = 'batch';
data.num  = [1 Inf];
data.help = {'Select the batch job files to be run.'};

opts.type = 'branch';
opts.name = 'Execute Batch Jobs';
opts.tag  = 'runbatch';
opts.val  = {data};
opts.prog = @runbatch;
opts.help = {[...
'This facility allows previously created batch jobs to be run. ',...
'These are simply created by the batch user interface ',...
'(which you are currently using).']};
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function runbatch(varargin)
jobs = varargin{1}.jobs;
for i=1:numel(jobs),
    spm_jobman('run',jobs{i});
end;
return;


function exec_jobman(arg1)
% A function to be compiled, which will run a batch job saved in job.xml
%_______________________________________________________________________
% John Ashburner $Id$ 
spm_defaults
if nargin==0, arg1 = 'job.xml'; end;
spm_jobman('run',arg1);

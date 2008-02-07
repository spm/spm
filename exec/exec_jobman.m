function exec_jobman(arg1)
% A function to be compiled, which will run a batch job saved in job.xml
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner 
% $Id: exec_jobman.m 1143 2008-02-07 19:33:33Z spm $ 

spm_defaults
path(path,spm('Dir'));
if nargin==0, arg1 = 'job.xml'; end;
spm_jobman('run',arg1);

function exec_jobman(arg1)
% A function to be compiled, which will run a batch job saved in job.xml
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner 
% $Id: exec_jobman.m 253 2005-10-13 15:31:34Z guillaume $ 

spm_defaults
path(path,spm('Dir'));
if nargin==0, arg1 = 'job.xml'; end;
spm_jobman('run',arg1);

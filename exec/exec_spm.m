function exec_spm(arg1)
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: exec_spm.m 253 2005-10-13 15:31:34Z guillaume $ 

path(path,spm('Dir'));
if nargin==0,
    spm;
else
    spm(arg1);
end;

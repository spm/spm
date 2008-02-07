function exec_spm(arg1)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: exec_spm.m 1143 2008-02-07 19:33:33Z spm $ 

path(path,spm('Dir'));
if nargin==0,
    spm;
else
    spm(arg1);
end;

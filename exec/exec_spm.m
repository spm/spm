function exec_spm(arg1)
%_______________________________________________________________________
% John Ashburner $Id$ 
path(path,spm('Dir'));
if nargin==0,
    spm;
else
    spm(arg1);
end;

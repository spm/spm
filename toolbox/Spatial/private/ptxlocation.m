function ptx = ptxlocation(nam)
% Location of a ptx file

pth = mfilename('fullpath');
[d,~,~] = fileparts(pth);
ptx = fullfile(d,"..","lib",parallel.gpu.ptxext);

if nargin>=1
    [~,nam,ext] = fileparts(nam);
    nam = [nam '.ptx'];
    ptx = fullfile(ptx,nam);
end


function M = off_read(filename)
% Read OFF-formatted data from disk
% FORMAT M = off_read(filename)
%
% filename - OFF-formatted file name
% M        - data structure
%__________________________________________________________________________
% 
% OFF Format Specification:
% https://en.wikipedia.org/wiki/OFF_(file_format)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: off_read.m 7676 2019-10-22 10:19:29Z guillaume $


fid = fopen(filename,'rt');
if fid == -1
    error('Cannot open %s.',filename);
end

M = struct('vertices',[],'faces',[]);

%-Read header
while true
    l = fgetl(fid);
    if isempty(l) || strcmp(l,'OFF') || l(1) == '#'
        continue;
    end
    n = sscanf(l,'%d %d %d');
    break;
end

%-Read data
M.vertices = fscanf(fid,'%f',[3 n(1)])';
M.faces    = fscanf(fid,'%d',[4 n(2)])';
M.faces    = M.faces(:,2:4) + 1;

fclose(fid);

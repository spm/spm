function M = ply_read(filename)
% Read PLY-formatted data from disk
% FORMAT M = ply_read(filename)
%
% filename - PLY-formatted file name
% M        - data structure
%__________________________________________________________________________
% 
% Stanford Triangle Format Specification:
% https://en.wikipedia.org/wiki/PLY_%28file_format%29
%__________________________________________________________________________
% Copyright (C) 2017-2019 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: ply_read.m 7676 2019-10-22 10:19:29Z guillaume $


fid = fopen(filename,'rt');
if fid == -1
    error('Cannot open %s.',filename);
end

M = struct('vertices',[],'faces',[]);

%-Read header
isvertex = false;
nval = 0;
while true
    l = fgetl(fid);
    if strcmp(l,'end_header'), break; end
    [l,r] = strtok(l);
    if strcmp(l,'element')
        [l,r] = strtok(r);
        switch l
            case 'vertex'
                nv = str2double(r);
                isvertex = true;
            case 'face'
                nf = str2double(r);
                isvertex = false;
        end
    elseif strcmp(l,'property')
        if isvertex
            nval = nval + 1;
        end
    end
end

%-Read data
M.vertices = fscanf(fid,'%f',[nval nv])';
M.vertices = M.vertices(:,1:3);
M.faces    = fscanf(fid,'%d',[4 nf])';
M.faces    = M.faces(:,2:4) + 1;

fclose(fid);

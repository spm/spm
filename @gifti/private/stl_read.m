function M = stl_read(filename)
% Read STL-formatted data from disk
% FORMAT M = stl_read(filename)
%
% filename - STL-formatted file name
% M        - data structure
%__________________________________________________________________________
% 
% STL Format Specification:
% https://en.wikipedia.org/wiki/STL_%28file_format%29
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: stl_read.m 7254 2018-02-05 17:51:37Z guillaume $


M = struct('vertices',[],'faces',[]);

%-Detect ASCII vs Binary STL
%--------------------------------------------------------------------------
fid = fopen(filename,'rt');
if fid == -1
    error('Cannot open %s.',filename);
end
hdr = fscanf(fid,'%c',80);
fclose(fid);

%-ASCII STL
%--------------------------------------------------------------------------
if strncmp(hdr,'solid',5)
    % solid name
    % facet normal ni nj nk
    %   outer loop
    %     vertex v1x v1y v1z
    %     vertex v2x v2y v2z
    %     vertex v3x v3y v3z
    %   endloop
    % endfacet
    % endsolid name
    fid = fopen(filename,'rt');
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        tline = strtrim(tline);
        switch strtok(tline)
            case {'solid','facet','outer','endloop','endfacet'}
            case 'vertex'
                tri = sscanf(tline,'vertex %f %f %f');
                M.vertices = [M.vertices tri];
        end
    end
    fclose(fid);
    M.vertices = M.vertices';
    M.faces = reshape(1:size(M.vertices,1),3,[])';
    return;
end

%-Binary STL
%--------------------------------------------------------------------------
fid = fopen(filename,'r','ieee-le');
% INT8[80] – Header
fread(fid,80,'uchar');
% UINT32 – Number of triangles
nf = fread(fid,1,'uint32');
M.vertices = zeros(3,3*nf);
M.faces = zeros(3,nf);
% foreach triangle
% REAL32[3] – Normal vector
% REAL32[3] – Vertex 1
% REAL32[3] – Vertex 2
% REAL32[3] – Vertex 3
% UINT16 – Attribute byte count
% end
for i=1:nf
    d = fread(fid,12,'float32');
    fread(fid,1,'uint16');
    M.vertices((1:9)+9*(i-1)) = d(4:end);
    M.faces((1:3)+3*(i-1)) = (1:3)+3*(i-1);
end
fclose(fid);
M.vertices = M.vertices';
M.faces = M.faces';

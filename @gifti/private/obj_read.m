function M = obj_read(filename)
% Read Wavefront OBJ-formatted data from disk
% FORMAT M = obj_read(filename)
%
% filename - OBJ-formatted file name
% M        - data structure
%__________________________________________________________________________
% 
% Wavefront OBJ Format Specification:
% https://en.wikipedia.org/wiki/Wavefront_.obj_file
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: obj_read.m 7001 2017-02-02 13:53:16Z guillaume $


fid = fopen(filename,'rt');
if fid == -1
    error('Cannot open %s.',filename);
end

M = struct('vertices',[],'faces',[]);

while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    if numel(l) < 1 || l(1) == '#', continue; end
    switch l(1)
        case 'v'
            switch l(2)
                case 't'
                case 'n'
                case 'p'
                otherwise
                    v = sscanf(l(2:end),'%f %f %f');
                    if numel(v) > 3, v = v(1:3); end
                    M.vertices(size(M.vertices,1)+1,:) = v;
            end
        case 'f'
            f = sscanf(l(2:end),'%d %d %d');
            if numel(f) ~= 3
                f = sscanf(l(2:end),'%d/%d %d/d %d/%d');
                if numel(f) ~= 6
                    f = sscanf(l(2:end),'%d//%d %d//d %d//%d');
                    if numel(f) ~= 6
                        f = sscanf(l(2:end),'%d/%d/%d %d/%d/%d %d/%d/%d');
                        if numel(f) == 9
                            f = f([1 4 7]);
                        else
                            fprintf('Not a triangle.\n');
                            continue;
                        end
                    else
                        f = f([1 3 5]);
                    end
                else
                    f = f([1 3 5]);
                end
            end
            i = find(f<0);
            if isempty(i), f(i) = size(M.vertices,1) + f(i); end
            M.faces(size(M.faces,1)+1,:) = f;
        case 'o'
            fprintf('Igoring named objects.\n');
        case 'g'
            fprintf('Igoring polygon groups.\n');
        case 's'
            fprintf('Igoring smooth shading.\n');
        otherwise
            if ~isempty(strmatch('mtllib',l)) || ~isempty(strmatch('usemtl',l))
                fprintf('Igoring materials.\n');
            else
                fprintf('Igoring line starting with %c.\n',l(1));
            end
    end
end

fclose(fid);

function M = mvtk_read(filename)
% Read VTK formatted data from disk
% FORMAT M = mvtk_read(filename)
%
% filename - VTK-formatted file name
% M        - data structure
%__________________________________________________________________________
% 
% VTK File Formats Specifications:
% http://www.vtk.org/VTK/img/file-formats.pdf
% 
% Requirements: zstream, base64decode
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: mvtk_read.m 7354 2018-06-22 10:44:22Z guillaume $


[pth,name,ext] = fileparts(filename);
switch ext
    case '.vtk'
        M = mvtk_read_legacy(filename);
    case {'.vti','.vtp','.vtr','.vts','.vtu'}
        % Serial vtkImageData (structured)
        % Serial vtkPolyData (unstructured)
        % Serial vtkRectilinearGrid (structured)
        % Serial vtkStructuredGrid (structured)
        % Serial vtkUnstructuredGrid (unstructured)
        M = mvtk_read_xml(filename);
    otherwise
        error('Unknown file format.');
end

%==========================================================================
% function M = mvtk_read_legacy(filename)
%==========================================================================
function M = mvtk_read_legacy(filename)

fid = fopen(filename,'rt');
if fid == -1
    error('Cannot open %s.',filename);
end

%- Part 1: file version and identifier
% # vtk DataFile Version 2.0
l = fgetl(fid);
if ~ischar(l), error('VTK file incomplete.'); end
if ~strncmpi(l,'# vtk DataFile Version',22)
    error('This is not a VTK formatted file.');
end

%- Part 2: header
l = fgetl(fid);
if ~ischar(l), error('VTK file incomplete.'); end

%- Part 3: file format
format = fgetl(fid);
if ~ismember(format,{'ASCII','BINARY'})
    error('Unknown file format.');
end

%- Part 4: dataset structure
data_attributes = false;
l = fgetl(fid);
if ~ischar(l), error('VTK file incomplete.'); end
[D,l] = strtok(l);
if ~strcmp(D,'DATASET'), error('Invalid VTK file.'); end
F = strtok(l(2:end));
switch F
    case 'STRUCTURED_POINTS'
        warning('Unsupported dataset format.');
        l = fgetl(fid);
        DIM = sscanf(l,'DIMENSIONS %d %d %d');
        l = fgetl(fid);
        ORIGIN = sscanf(l,'ORIGIN %f %f %f');
        l = fgetl(fid);
        SPACING = sscanf(l,'SPACING %f %f %f');
    case 'STRUCTURED_GRID'
        warning('Unsupported dataset format.');
        l = fgetl(fid);
        DIM = sscanf(l,'DIMENSIONS %d %d %d');
        l = fgetl(fid);
        % assume float and n = prod(DIM)
        PTS = textscan(fid,'%f %f %f\n',prod(DIM),'CollectOutput',true);
        PTS = PTS{1};
    case 'RECTILINEAR_GRID'
        warning('Unsupported dataset format.');
        l = fgetl(fid);
        DIM = sscanf(l,'DIMENSIONS %d %d %d');
        l = fgetl(fid);
        XCOORDS = textscan(fid,'%f',DIM(1),'CollectOutput',true);
        XCOORDS = XCOORDS{1};
        l = fgetl(fid);
        YCOORDS = textscan(fid,'%f',DIM(2),'CollectOutput',true);
        YCOORDS = YCOORDS{1};
        l = fgetl(fid);
        ZCOORDS = textscan(fid,'%f',DIM(3),'CollectOutput',true);
        ZCOORDS = ZCOORDS{1};
    case 'POLYDATA'
        while true
            l = fgetl(fid);
            if ~ischar(l), break; end
            [D,l] = strtok(l);
            switch D
                case 'POINTS'
                    [N,l] = strtok(l(2:end)); % l still contains dataType
                    N = str2double(N);
                    M.vertices = textscan(fid,'%f %f %f\n',N,'CollectOutput',true);
                    M.vertices = M.vertices{1};
                case 'POLYGONS'
                    [N,l] = strtok(l(2:end));
                    N = str2double(N);
                    S = strtok(l);
                    S = str2double(S);
                    if 4*N ~= S, error('Unsupported dataset format.'); end
                    M.faces = textscan(fid,'3 %d %d %d\n',N,'CollectOutput',true);
                    M.faces = M.faces{1} + 1;
                case {'VERTICES','LINES','TRIANGLE_STRIPS'}
                    error('Unsupported data type.');
                case {'POINT_DATA','CELL_DATA'}
                    data_attributes = true;
                    [N,l] = strtok(l(2:end));
                    N = str2double(N);
                    break;
                case ''
                otherwise
                    error('Invalid VTK file.');
            end
        end
    case {'UNSTRUCTURED_GRID','FIELD'}
        error('Unsupported data type.');
    otherwise
        error('Invalid VTK file.');
end

%- Part 5: dataset attributes (POINT_DATA and CELL_DATA)
if data_attributes
    %l = fgetl(fid); % {POINT_DATA,CELL_DATA} N
    l = fgetl(fid); % SCALARS dataName dataType numComp
    [P,l] = strtok(l);
    if strcmp(P,'SCALARS')
        [S,l] = strtok(l(2:end));
        [S,l] = strtok(l(2:end));
        S = strtok(l(2:end)); S = str2double(S);
        l = fgetl(fid); % LOOKUP_TABLE default
        fmt = repmat('%f ',1,S);
        fmt = [fmt(1:end-1) '\n'];
        M.cdata = textscan(fid,fmt,N,'CollectOutput',true);
        M.cdata = M.cdata{1};
    end
end

fclose(fid);

%==========================================================================
% function M = mvtk_read_xml(filename)
%==========================================================================
function M = mvtk_read_xml(filename)

try
    X = xmltree(filename);
catch
    error('Cannot parse file %s.',filename);
end

warning('Unsupported file format.');
M = struct([]);

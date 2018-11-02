function this = read_freesurfer_file(filename)
% Low level reader of FreeSurfer file
% FORMAT this = read_freesurfer_file(filename)
% filename    - FreeSurfer file
%
% Read ASCII triangle surface file and part of binary mgh file.
% See https://surfer.nmr.mgh.harvard.edu/fswiki/FileFormats
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: read_freesurfer_file.m 7471 2018-11-02 11:14:39Z guillaume $


[p,n,e] = fileparts(filename);
switch lower(e)
    case {'.asc','.srf'}
        this = read_fs_file_ascii(filename);
    case {'.mgh','.mgz'}
        this.cdata = read_fs_file_mgh(filename);
    case {'.pial','.white','.inflated','.nofix','.orig','.smoothwm',...
            '.sphere','.reg','.surf'}
        this = read_fs_file_mesh(filename);
    case {'.curv','.area','.sulc'}
        this = read_fs_file_cdata(filename);
    otherwise
        error('Unknown file format.');
end


function this = read_fs_file_ascii(filename)
fid = fopen(filename,'rt');
if fid == -1, error('Cannot open "%s".',filename); end

fgetl(fid); % #!ascii 
N = fscanf(fid,'%d',2);
this.vertices = fscanf(fid,'%f %f %f %d',[4 N(1)])';
this.faces    = fscanf(fid,'%d %d %d %d',[4 N(2)])';

fclose(fid);

this.vertices = this.vertices(:,1:3);
this.faces    = this.faces(:,1:3) + 1;


function [dat,hdr] = read_fs_file_mgh(filename)
[p,n,e] = fileparts(filename);
if strcmpi(e,'.mgz')
    filename = char(gunzip(filename,tempname));
    c = onCleanup(@()rmdir(fileparts(filename),'s'));
end

fid = fopen(filename,'r','ieee-be');
if fid == -1, error('Cannot open "%s".',filename); end

% version integer current value is 1
hdr.version = fread(fid,1,'int32');
% width integer first dimension of the image buffer (fastest)
% height integer second dimension of the image buffer (2nd fastest)
% depth integer slowest dimension when reading the buffer
% nframes integer number of scalar components contained in the buffer
hdr.dim = fread(fid,4,'int32')';
% type integer data type of the image buffer: 0 (UCHAR), 4 (SHORT), 1 (INT) or 3 (FLOAT)
hdr.type = fread(fid,1,'int32');
% dof integer degrees of freedom (where appropriate)
hdr.dof = fread(fid,1,'int32');
% goodRASFlag short if true, the direction cosines are in the header
% if false, a CORONAL orientation will be assumed
hdr.goodRASFlag = fread(fid,1,'int16');
hdr.spacing = [1 1 1]';
hdr.xyz = [-1 0 0;0 0 1;0 -1 0];
hdr.c = [0 0 0]';
if hdr.goodRASFlag
    % spacing X float spacing in the X direction (ranging [0...width]) - default is 1
    % spacing Y float spacing in the Y direction (ranging [0...height]) - default is 1
    % spacing Z float spacing in the Z direction (ranging [0...depth]) - default is 1
    hdr.spacing = fread(fid,3,'float32');
    % xr float default is -1
    % xa float default is 0
    % xs float default is 0
    % yr float default is 0
    % ya float default is 0
    % ys float default is -1
    % zr float default is 0
    % za float default is 1
    % zs float default is 0
    hdr.xyz = reshape(fread(fid,9,'float32'),[3 3]);
    % cr float default is 0
    % ca float default is 0
    % cs float default is 0
    hdr.c = fread(fid,3,'float32');
end
% The image data starts at a specified offset from the beginning of the file,
% which is currently byte 284 (bytes 0-283 are header)
fseek(fid,284,'bof');
dt = {'uchar', 'int32', NaN, 'float32', 'int16'};
dat = fread(fid,prod(hdr.dim),dt{hdr.type+1});
dat = reshape(dat,hdr.dim);

% Immediately after the data buffer can be found optional data structures
% TR float milliseconds
% FlipAngle float radians
% TE float milliseconds
% TI float milliseconds
if ~feof(fid)
    hdr.params = fread(fid,4,'float32');
end
% tags variable length char strings
if ~feof(fid)
    % hdr.tags = fread(fid,Inf,'uchar');
end

fclose(fid);


function this = read_fs_file_mesh(filename)
fid = fopen(filename,'r','ieee-be');
if fid == -1, error('Cannot open "%s".',filename); end

magic = fread(fid,3,'uchar');
magic = [65536 256 1]*magic;
if magic ~= 16777214
    fclose(fid);
    error('File format not supported.');
end
fgets(fid); fgets(fid);
nv = fread(fid,1,'int32');
nf = fread(fid,1,'int32');
this.vertices = fread(fid,3*nv,'float32');
this.vertices = reshape(this.vertices,3,nv)';
this.faces = fread(fid,3*nf,'int32');
this.faces = reshape(this.faces,3,nf)' + 1;

fclose(fid);


function this = read_fs_file_cdata(filename)
fid = fopen(filename,'r','ieee-be');
if fid == -1, error('Cannot open "%s".',filename); end

magic = fread(fid,3,'uchar');
magic = [65536 256 1]*magic;
if magic ~= 16777215
    fclose(fid);
    error('File format not supported.');
end

nv = fread(fid,1,'int32');
nf = fread(fid,1,'int32');
n  = fread(fid,1,'int32');
this.cdata = fread(fid,[n nv],'float32')';

fclose(fid) ;

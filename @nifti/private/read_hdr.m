function vol = read_hdr(fname)
% Get a variety of information from a NIFTI-1 header.
% FORMAT vol = read_hdr(fname)
% fname - filename of image
% vol   - various bits of information
% _______________________________________________________________________
% $Id$

persistent dict
if isempty(dict), dict = getdict; end;

[pth,nam,ext] = fileparts(fname);
switch ext
case '.hdr'
    hname = fullfile(pth,[nam '.hdr']);
case '.HDR'
    hname = fullfile(pth,[nam '.HDR']);
case '.img'
    hname = fullfile(pth,[nam '.hdr']);
case '.IMG'
    hname = fullfile(pth,[nam '.HDR']);
case '.nii'
    hname = fullfile(pth,[nam '.nii']);
case '.NII'
    hname = fullfile(pth,[nam '.NII']);
otherwise
    hname = fullfile(pth,[nam ext]);
end
[hdr,be] = read_hdr_raw(hname);
hdr      = mayo2nifti1(hdr);
d  = getdict;
dt = [];
for i=1:length(d.dtype)
    if hdr.datatype == d.dtype(i).code
        dt = d.dtype(i);
        break;
    end;
end;
if isempty(dt)
    error(['Unrecognised datatype (' num2str(double(hdr.datatype)) ') for "' fname '.'] );
end
if isfield(hdr,'magic')
    switch deblank(hdr.magic)
    case {'n+1'}
        iname = hname;
        if hdr.vox_offset < hdr.sizeof_hdr
            error(['Bad vox_offset (' num2str(double(hdr.vox_offset)) ') for "' fname '.'] );
        end
    case {'ni1'}
        if strcmp(ext,lower(ext)),
            iname = fullfile(pth,[nam '.img']);
        else,
            iname = fullfile(pth,[nam '.IMG']);
        end;
    otherwise
        error(['Bad magic (' hdr.magic ') for "' fname '.'] );
    end
else
    if strcmp(ext,lower(ext)),
        iname = fullfile(pth,[nam '.img']);
    else,
        iname = fullfile(pth,[nam '.IMG']);
    end;
end
if rem(double(hdr.vox_offset),dt.size)
   error(['Bad alignment of voxels (' num2str(double(hdr.vox_offset)) '/' num2str(double(dt.size)) ') for "' fname '.'] );
end;

vol      = struct('hdr',hdr,'be',be,'hname',hname,'iname',iname);
return

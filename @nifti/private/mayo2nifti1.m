function hdr = mayo2nifti1(ohdr)
% Convert from an ANALYZE to a NIFTI-1 header
% _______________________________________________________________________
% $Id$

if isfield(ohdr,'magic'),
    hdr = ohdr;
    return;
end;
hdr            = empty_hdr;
hdr.dim        = ohdr.dim;
hdr.datatype   = ohdr.datatype;
hdr.bitpix     = ohdr.bitpix;
hdr.pixdim     = ohdr.pixdim;
hdr.vox_offset = ohdr.vox_offset;
hdr.scl_slope  = ohdr.roi_scale;
hdr.scl_inter  = ohdr.funused1;
hdr.descrip    = ohdr.descrip;
hdr.aux_file   = ohdr.aux_file;
hdr.glmax      = ohdr.glmax;
hdr.glmin      = ohdr.glmin;
hdr.cal_max    = ohdr.cal_max;
hdr.cal_min    = ohdr.cal_min;
hdr.magic      = 'ni1';

%mat            = decode_qform0(ohdr);
%Make this more SPM format specific...
if any(ohdr.origin(1:3)), origin = double(ohdr.origin(1:3));
else,                     origin = (double(ohdr.dim(2:4))+1)/2; end;
vox    = double(ohdr.pixdim(2:4));
if all(vox == 0), vox = [1 1 1]; end;
off    = -vox.*origin;
mat    = [vox(1) 0 0 off(1) ; 0 vox(2) 0 off(2) ; 0 0 vox(3) off(3) ; 0 0 0 1];
if spm_flip_analyze_images, mat = diag([-1 1 1 1])*mat; end;

hdr            = encode_qform0(mat,hdr);
mat            = mat*[eye(4,3) [1 1 1 1]'];
hdr.srow_x     = mat(1,:);
hdr.srow_y     = mat(2,:);
hdr.srow_z     = mat(3,:);
return;

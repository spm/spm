function [N,cdf] = spm_mnc2nifti(fname,opts)
% Import MINC images into NIfTI
% FORMAT spm_mnc2nifti(fname)
% fname - a MINC filename
% opts  - options structure
%
% N     - NIfTI object (written in current directory)
% cdf   - NetCDF data structure
% 
% The MINC file format was developed by Peter Neelin at the Montreal
% Neurological Institute, and is based upon the NetCDF libraries.
% The NetCDF documentation specifically recommends that people do not
% write their own libraries for accessing the data.  This suggestion
% was ignored.
%
% See: http://en.wikibooks.org/wiki/MINC
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_mnc2nifti.m 4927 2012-09-14 16:15:10Z ged $


if nargin==1
    opts = struct('dtype',4, 'ext',spm_file_ext);
else
    if opts.ext(1) ~= '.', opts.ext = ['.' opts.ext]; end
end

cdf = spm_read_netcdf(fname);
if isempty(cdf)
    error(['"' fname '" does not appear to be MINC.']);
end

%-Create file_array object
%--------------------------------------------------------------------------
idat = file_array;

d_types     = [2   2   512   768     16   64 ; 256  256  4      8      16   64];
%dsizes     = [1   1   2     4       4    8 ];
mxs         = [255 255 65535 2^32-1  Inf  Inf; 127  127  32767  2^31-1 Inf  Inf];
mns         = [0   0   0     0      -Inf -Inf;-128 -128 -32768 -2^31  -Inf -Inf];
%space_names= {'xspace','yspace','zspace'};
img         = findvar(cdf.var_array,'image');
nd          = length(img.dimid);
idat.fname  = fname;
idat.dim    = fliplr(cat(2,cdf.dim_array(:).dim_length));
signed      = findvar(img.vatt_array,'signtype');
signed      = strcmp(signed.val,'signed__');
idat.dtype  = [d_types(signed+1,img.nc_type) 1];
range       = [mns(signed+1,img.nc_type) mxs(signed+1,img.nc_type)];
idat.offset = img.begin;


if img.nc_type <=4
    tmp = findvar(img.vatt_array,'valid_range');
    if isempty(tmp)
        disp(['Can''t get valid_range for "' fname '" - having to guess']);
    else
        range = tmp.val;
    end

    fp   = fopen(fname,'r','ieee-be');
    imax = get_imax(fp, cdf, 'image-max', 1, idat.dim);
    imin = get_imax(fp, cdf, 'image-min', 0, idat.dim);
    fclose(fp);

    scale = (imax-imin)/(range(2)-range(1));
    dcoff = imin-range(1)*scale;
else
    scale = 1;
    dcoff = 0;
end

%-Extract affine transformation from voxel to world co-ordinates
%--------------------------------------------------------------------------
step  = [1 1 1];
start = [0 0 0]';
dircos = eye(3);
for j=1:3
    nam    = cdf.dim_array(img.dimid(nd+1-j)).name;
    space  = findvar(cdf.var_array,nam);
    tmp    = findvar(space.vatt_array,'step');
    if ~isempty(tmp), step(j) = tmp.val; end
    tmp    = findvar(space.vatt_array,'start');
    if ~isempty(tmp), start(j) = tmp.val; else  start(j) = -dim(j)/2*step(j); end
    tmp    = findvar(space.vatt_array,'direction_cosines');
    if ~isempty(tmp)
        if tmp.nc_type == 6
            dircos(:,j) = tmp.val(:);
        elseif tmp.nc_type == 2
            dircos(:,j) = sscanf(tmp.val,'%g');
        end
    end
end
shiftm = [1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1];
mat    = [[dircos*diag(step) dircos*start] ; [0 0 0 1]] * shiftm;

%-Create NIfTI object
%--------------------------------------------------------------------------
N          = nifti;
dat        = file_array;
nam        = spm_file(fname,'basename');
dat.fname  = fullfile(pwd,[nam opts.ext]);
dat.dim    = idat.dim;
dat.dtype  = [opts.dtype spm_platform('bigend')];
dat.offset = 0;

if ~spm_type(opts.dtype,'intt')
    dat.scl_slope = 1;
    dat.scl_inter = 0;
else
    mn = Inf;
    mx = -Inf;
    for i6=1:size(idat,6)
        for i5=1:size(idat,5)
            for i4=1:size(idat,4)
                for i3=1:size(idat,3)
                    if size(scale,3)==1
                        scale1 = scale(:,:,1,i4,i5,i6);
                        dcoff1 = dcoff(:,:,1,i4,i5,i6);
                    else
                        scale1 = scale(:,:,i3,i4,i5,i6);
                        dcoff1 = dcoff(:,:,i3,i4,i5,i6);
                    end
                    if numel(scale1)==1
                        img = double(idat(:,:,i3,i4,i5,i6))*scale1 + dcoff1;
                    elseif size(scale1,1)>1 && size(scale1,2)>1
                        img = double(idat(:,:,i3,i4,i5,i6)).*scale1 + dcoff1;
                    elseif size(scale1,1)==1
                        img = double(idat(:,:,i3,i4,i5,i6)).*repmat(scale1,[size(idat,1) 1]) +...
                                                             repmat(dcoff1,[size(idat,1) 1]);
                    else
                        img = double(idat(:,:,i3,i4,i5,i6)).*repmat(scale1,[1 size(idat,2)]) +...
                                                             repmat(dcoff1,[1 size(idat,2)]);
                    end
                    img = img(isfinite(img));
                    if ~isempty(img)
                        mx  = max(mx,max(img));
                        mn  = min(mn,min(img));
                    end
                end
            end
        end
    end
    dat.scl_slope = mx/spm_type(opts.dtype,'maxval');
    if spm_type(opts.dtype,'minval')~=0
        dat.scl_slope = max(dat.scl_slope,mn/spm_type(opts.dtype,'minval'));
    end
    dat.scl_inter = 0;
end

N.dat   = dat;
flp     = false;
if det(mat)>0
    flp = true;
    mat = mat*[diag([-1 1 1]) [size(dat,1)+1 0 0]' ; 0 0 0 1];
end
N.mat  = mat;
N.mat0 = mat;
create(N);
for i6=1:size(idat,6)
    for i5=1:size(idat,5)
        for i4=1:size(idat,4)
            for i3=1:size(idat,3)
                if size(scale,3)==1
                    scale1 = scale(:,:,1,i4,i5,i6);
                    dcoff1 = dcoff(:,:,1,i5,i5,i6);
                else
                    scale1 = scale(:,:,i3,i4,i5,i6);
                    dcoff1 = dcoff(:,:,i3,i4,i5,i6);
                end
                if numel(scale1)==1
                    slice = double(idat(:,:,i3,i4,i5,i6))*scale1 + dcoff1;
                elseif size(scale1,1)>1 && size(scale1,2)>1
                    slice = double(idat(:,:,i3,i4,i5,i6)).*scale1 + dcoff1;
                elseif size(scale1,1)==1
                    slice = double(idat(:,:,i3,i4,i5,i6)).*repmat(scale1,[size(idat,1) 1]) +...
                                                           repmat(dcoff1,[size(idat,1) 1]);
                else
                    slice = double(idat(:,:,i3,i4,i5,i6)).*repmat(scale1,[1 size(idat,2)]) +...
                                                           repmat(dcoff1,[1 size(idat,2)]);
                end
                if flp, slice = flipud(slice); end
                N.dat(:,:,i3,i4,i5,i6) = slice;
            end
        end
    end
end


%==========================================================================
% function var = findvar(varlist, name)
%==========================================================================
function var = findvar(varlist, name)
% Finds the structure in a list of structures that has a name element
% matching the second argument.
for i=1:numel(varlist)
    if strcmp(varlist(i).name,name)
        var = varlist(i);
        return;
    end
end
var = [];
%error(['Can''t find "' name '".']);


%==========================================================================
% function str = dtypestr(i)
%==========================================================================
function str = dtypestr(i)
% Returns a string appropriate for reading or writing the CDF data-type.
types = char('uint8','uint8','int16','int32','float32','float64');
str   = deblank(types(i,:));


%==========================================================================
% function imax = get_imax(fp, cdf, strng, def, dim)
%==========================================================================
function imax = get_imax(fp, cdf, strng, def, dim)
dim  = fliplr(dim);
str  = findvar(cdf.var_array,strng);

if ~isempty(str) && str.nc_type == 6
    fseek(fp,str.begin,'bof');
    nel  = str.vsize/(spm_type(dtypestr(str.nc_type),'bits')/8);
    imax = fread(fp,nel,dtypestr(str.nc_type))';
    resh = ones(1,numel(dim));
    resh(numel(dim)+1-str.dimid) = dim(str.dimid);
    imax = reshape(imax,resh);
else
    imax = def;
end

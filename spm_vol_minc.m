function V=spm_vol_minc(fname)
% Get header information etc. for MINC images.
%  FORMAT V = spm_vol(P)
%  P - a MINC filename.
%  V - a structure containing image volume information.
%
% The elements of V are described by the help for spm_vol, except for
% an additional field (V.cdf) that contains the NetCDF header
% information.
%
% The MINC file format was developed by Peter Neelin at the Montréal
% Neurological Institute, and is based upon the NetCDF libraries.
% The NetCDF documentation specifically recommends that people do not
% write their own libraries for accessing the data.  However, this
% was necessary in this case, in order to attempt to get MINC images into
% a format that the spm image reading utilities could use.
% _______________________________________________________________________
%  %W% John Ashburner %E%

cdf = read_netcdf(fname);
if isempty(cdf), V=[]; return; end;


d_types = [2 2 4 8 16 64];
dsizes  = [1 1 2 4  4  8];

img   = findvar(cdf.var_array,'image');

% Extract the dimensions.  This is stuff I don't fully understand
% so there may be problems.
%-----------------------------------------------------------------------
for i=1:prod(size(img.dimid)),
	dim0(i) = cdf.dim_array(img.dimid(i));
end;

dim = fliplr(cat(2,dim0.dim_length));
dim = dim(1:min(size(dim,2),3));
dim = [dim ones(1,3-size(dim,2))];

datatype = d_types(img.nc_type);
signed   = findvar(img.vatt_array,'signtype');
signed   = strcmp(signed.val,'signed__');
is_flt   = 0;
switch datatype,
	case {2},
		if signed,
			datatype = datatype + 128;
			range = [-2^7 2^7-1]';
		else,
			range = [0 2^8-1]';
		end;
	case {4},
		if signed,
			range = [-2^15 2^15-1]';
		else,
			datatype = datatype + 128;
			range = [0 2^16-1]';
		end;
	case {8},
		if signed,
			range = [-2^31 2^31-1]';
		else,
			datatype = datatype + 128;
			range = [0 2^32-1]';
		end;
	otherwise,
		is_flt = 1;
end;
if ~spm_platform('bigend') & datatype~=2, datatype = datatype*256; end;

dim   = [dim(1:3) datatype];

if ~is_flt,
	tmp = findvar(img.vatt_array,'valid_range');
	if isempty(tmp),
		disp(['Can''t get valid_range for "' fname '" - having to guess']);
	else
		range = tmp.val;
	end;

	fp    = fopen(fname,'r','ieee-be');
	imax  = findvar(img.vatt_array,'image-max');
	if ~isempty(imax) & imax.nc_type == 2,
		imax = findvar(cdf.var_array,imax.val(5:end));
		fseek(fp,imax.begin,'bof');
		ddim = dim(3);
		imax = fread(fp,prod(ddim),dtypestr(imax.nc_type))';
		imax = reshape(imax,[ddim 1 1]);
		imax = imax(:,1,1)';
	else,
		imax = ones(1,dim(3));
		disp(['Can''t get imax for "' fname '" - guessing it is one']);
	end;
	imin  = findvar(img.vatt_array,'image-min');
	if ~isempty(imin) & imin.nc_type == 2,
		imin = findvar(cdf.var_array,imin.val(5:end));
		fseek(fp,imin.begin,'bof');
		ddim = dim(3);
		imin = fread(fp,prod(ddim),dtypestr(imin.nc_type))';
		imin = reshape(imin,[ddim 1 1]);
		imin = imin(:,1,1)';
	else,
		imin = zeros(1,dim(3));
		disp(['Can''t get imax for "' fname '" - guessing it is zero']);
	end;
	fclose(fp);
	scale = (imax-imin)/(range(2)-range(1));
	dcoff = imin-range(1)*scale;
else,
	scale = ones(1,dim(3));
	dcoff = zeros(1,dim(3));
end;

off   = img.begin;
psize = dsizes(img.nc_type)*prod(dim(1:2));
off   = cumsum(ones(1,dim(3))*psize)-psize+off;

pinfo = [scale ; dcoff ; off];

% Extract affine transformation from voxel to world co-ordinates
%-----------------------------------------------------------------------
step  = [1 1 1];
start = [0 0 0]';
dircos = eye(3);

for j=1:3,
	nam    = cdf.dim_array(4-j).name;
	space  = findvar(cdf.var_array,nam);
	tmp    = findvar(space.vatt_array,'step');
	if ~isempty(tmp), step(j) = tmp.val; end;
	tmp    = findvar(space.vatt_array,'start');
	if ~isempty(tmp), start(j) = tmp.val; else, start(j) = -dim(j)/2*step(j); end;
	tmp    = findvar(space.vatt_array,'direction_cosines');
	if ~isempty(tmp), dircos(:,j) = tmp.val; end;
end;

shiftm = [1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1];
mat    = [[dircos*diag(step) dircos*start] ; [0 0 0 1]] * shiftm;

% Because there are not yet any routines to write the matrix information
% to MINC files, any changes to the matrix values will be made to `.mat'
% files.  The values in the `.mat' files should override the values from the
% headers.
matname = [spm_str_manip(fname,'sd') '.mat'];
if (exist(matname) == 2)
	load(matname,'M');
	if (exist('M') == 1)
		mat = M;
	end
end

V      = struct('fname',fname,'dim',dim,'mat',mat,'pinfo',pinfo,'cdf',cdf);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function var = findvar(varlist, name)
% Finds the structure in a list of structures that has a name element
% matching the second argument.
for i=1:prod(size(varlist)),
	if strcmp(varlist(i).name,name),
		var = varlist(i);
		return;
	end;
end;
var = [];
%error(['Can''t find "' name '".']);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function cdf = read_netcdf(fname)
% Read the header information from a CDF file into a data structure.
dsiz     = [1 1 2 4 4 8];
fp=fopen(fname,'r','ieee-be');
if fp==-1,
	cdf = [];
	return;
end;

% Return null if not a CDF file.
%-----------------------------------------------------------------------
mgc = fread(fp,4,'uchar')';
if ~all(['CDF' 1] == mgc),
	cdf = [];
	fclose(fp);
	return;
end

% I've no idea what this is for
numrecs = fread(fp,1,'uint32');

cdf = struct('numrecs',numrecs,'dim_array',[], 'gatt_array',[], 'var_array', []);

dt = fread(fp,1,'uint32');
if dt == 10,
	% Dimensions
	nelem = fread(fp,1,'uint32');
	for j=1:nelem,
		str   = readname(fp);
		dim_length = fread(fp,1,'uint32');
		cdf.dim_array(j).name       = str;
		cdf.dim_array(j).dim_length = dim_length;
	end;
	dt   = fread(fp,1,'uint32');
end

if dt == 12,
	% Attributes
	nelem = fread(fp,1,'uint32');
	for j=1:nelem,
		str    = readname(fp);
		nc_type= fread(fp,1,'uint32');
		nnelem = fread(fp,1,'uint32');
		val    = fread(fp,nnelem,dtypestr(nc_type));
		if nc_type == 2, val = deblank([val' ' ']); end
		padding= fread(fp,ceil(nnelem*dsiz(nc_type)/4)*4-nnelem*dsiz(nc_type),'uchar');
		cdf.gatt_array(j).name    = str;
		cdf.gatt_array(j).nc_type = nc_type;
		cdf.gatt_array(j).val     = val;
	end;
	dt   = fread(fp,1,'uint32');
end

if dt == 11,
	% Variables
	nelem = fread(fp,1,'uint32');
	for j=1:nelem,
		str    = readname(fp);
		nnelem = fread(fp,1,'uint32');
		val    = fread(fp,nnelem,'uint32');
		cdf.var_array(j).name    = str;
		cdf.var_array(j).dimid   = val+1;
		cdf.var_array(j).nc_type = 0;
		cdf.var_array(j).vsize   = 0;
		cdf.var_array(j).begin   = 0;
		dt0    = fread(fp,1,'uint32');
		if dt0 == 12,
			nelem0 = fread(fp,1,'uint32');
			for jj=1:nelem0,
				str    = readname(fp);
				nc_type= fread(fp,1,'uint32');
				nnelem = fread(fp,1,'uint32');
				val    = fread(fp,nnelem,dtypestr(nc_type));
				if nc_type == 2, val = deblank([val' ' ']); end
				padding= fread(fp,...
					ceil(nnelem*dsiz(nc_type)/4)*4-nnelem*dsiz(nc_type),'uchar');
				cdf.var_array(j).vatt_array(jj).name    = str;
				cdf.var_array(j).vatt_array(jj).nc_type = nc_type;
				cdf.var_array(j).vatt_array(jj).val     = val;
			end;
			dt0    = fread(fp,1,'uint32');
		end;
		cdf.var_array(j).nc_type  = dt0;
		cdf.var_array(j).vsize = fread(fp,1,'uint32');
		cdf.var_array(j).begin = fread(fp,1,'uint32');
	end;
	dt   = fread(fp,1,'uint32');
end;

fclose(fp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function str = dtypestr(i)
% Returns a string appropriate for reading or writing the CDF data-type.
types = str2mat('uint8','uint8','int16','int32','float','double');
str   = deblank(types(i,:));
return;
%_______________________________________________________________________

%_______________________________________________________________________
function name = readname(fp)
% Extracts a name from a CDF file pointed to at the right location by
% fp.
stlen  = fread(fp,1,'uint32');
name   = deblank([fread(fp,stlen,'uchar')' ' ']);
padding= fread(fp,ceil(stlen/4)*4-stlen,'uchar');
return;
%_______________________________________________________________________

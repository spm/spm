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
if ~strcmp(dim0(end).name,'xspace'),
	error('Fastest dimension does not appear to be X');
end;
if ~strcmp(dim0(end-1).name,'yspace'),
	error('Second fastest dimension does not appear to be Y');
end;
dim = [dim0(end).dim_length dim0(end-1).dim_length];
if ~strcmp(dim0(end-2).name,'zspace'),
	dim(3) = 1;
else
	dim(3) = dim0(end-2).dim_length;
end;


dim   = [dim(1:3) d_types(img.nc_type)];
if ~bigend & dim(4)~=2, dim(4) = dim(4)*256; end;

% This is something else I don't understand.  I thought that the
% stuff that is now commented out should work.  It does for most
% cases - except for some t-statistic images I tried.
%-----------------------------------------------------------------------
range = findvar(img.vatt_array,'valid_range');
range = range.val;

fp    = fopen(fname,'r','ieee-be');
imax  = findvar(img.vatt_array,'image-max');
if imax.nc_type == 2,
	imax = findvar(cdf.var_array,imax.val(5:end));
	fseek(fp,imax.begin,'bof');
	% For some reason this does not work....
	%if ~strcmp(dim0(imax.dimid(end)).name,'zspace')
	%	fclose(fp);
	%	error('Problem with dimensions');
	%end;
	%ddim = fliplr(cat(1,dim0(imax.dimid).dim_length));
	ddim = dim(3);
	imax = fread(fp,prod(ddim),dtypestr(imax.nc_type))';
	imax = reshape(imax,[ddim 1 1]);
	imax = imax(:,1,1)';
end;
imin  = findvar(img.vatt_array,'image-min');
if imin.nc_type == 2,
	imin = findvar(cdf.var_array,imin.val(5:end));
	fseek(fp,imin.begin,'bof');
	%if ~strcmp(dim0(imin.dimid(end)).name,'zspace')
	%	fclose(fp);
	%	V = [];
	%	return;
	%end;
	%ddim = fliplr(cat(1,dim0(imin.dimid).dim_length));
	ddim = dim(3);
	imin = fread(fp,prod(ddim),dtypestr(imin.nc_type))';
	imin = reshape(imin,[ddim 1 1]);
	imin = imin(:,1,1)';
end;
fclose(fp);

off   = img.begin;
psize = dsizes(img.nc_type)*prod(dim(1:2));
off   = cumsum(ones(1,dim(3))*psize)-psize+off;

scale = (imax-imin)/(range(2)-range(1));
dcoff = imin-range(1)*scale;
pinfo = [scale ; dcoff ; off];

% Extract affine transformation from voxel to world co-ordinates
%-----------------------------------------------------------------------
step  = [0 0 0];
start = [0 0 0]';
dircos = eye(3);
space  = findvar(cdf.var_array,'zspace');
tmp    = findvar(space.vatt_array,'step');
if ~isempty(tmp), step(3) = tmp.val; else, step(3) = 1; end;
tmp    = findvar(space.vatt_array,'start');
if ~isempty(tmp), start(3) = tmp.val; else, start(3) = -dim(3)/2*step(3); end;
tmp    = findvar(space.vatt_array,'direction_cosines');
if ~isempty(tmp), dircos(:,3) = tmp.val; end;

space  = findvar(cdf.var_array,'yspace');
tmp    = findvar(space.vatt_array,'step');
if ~isempty(tmp), step(2) = tmp.val; else, step(2) = 1; end;
tmp    = findvar(space.vatt_array,'start');
if ~isempty(tmp), start(2) = tmp.val; else, start(2) = -dim(2)/2*step(2); end;
tmp    = findvar(space.vatt_array,'direction_cosines');
if ~isempty(tmp), dircos(:,2) = tmp.val; end;

space  = findvar(cdf.var_array,'xspace');
tmp    = findvar(space.vatt_array,'step');
if ~isempty(tmp), step(1) = tmp.val; else, step(1) = 1; end;
tmp    = findvar(space.vatt_array,'start');
if ~isempty(tmp), start(1) = tmp.val; else, start(1) = -dim(1)/2*step(1); end;
tmp    = findvar(space.vatt_array,'direction_cosines');
if ~isempty(tmp), dircos(:,1) = tmp.val; end;

mat    = [[dircos*diag(step) dircos*start] ; [0 0 0 1]];

% Because there are not yet any routines to write the matrix information
% to MINC files, any changes to the matrix values will be made to `.mat'
% files.  The values in the `.mat' files should override the values from the
% headers.
matname = [spm_str_manip(fname,'sd') '.mat'];
if (exist(matname) == 2)
	load(matname);
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
		cdf.dim_array(j) = struct('name',str, 'dim_length', dim_length);
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
		cdf.gatt_array(j) = struct('name',str, 'nc_type',nc_type, 'val',val);
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
		cdf.var_array(j) = struct('name',str, 'dimid',val+1,...
			'vatt_array',[], 'nc_type',0, 'vsize', 0, 'begin', 0);
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
				cdf.var_array(j).vatt_array(jj) = struct(...
					'name',str, 'nc_type',nc_type, 'val',val);
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

%_______________________________________________________________________
function bend = bigend
% Checks to see if the computer is big-endian.  I don't know about some
% of the architectures - so it may be worth checking the code.
computers = str2mat('PCWIN','MAC2','SUN4','SOL2','HP700','SGI','SGI64','IBM_RS',...
			'ALPHA','AXP_VMSG','AXP_VMSIEEE','LNX86','VAX_VMSG','VAX_VMSD');
endians = [NaN NaN 1 1 NaN 0 0 NaN 0 Inf 0 NaN Inf Inf];
c=computer;
bend = NaN;
for i=1:size(computers,1),
	if strcmp(c,deblank(computers(i,:))),
		bend = endians(i);
		break;
	end;
end;
if ~finite(bend),
	if isnan(bend),
		error(['I don''t know if "' c '" is big-endian.']);
	else,
		error(['I don''t think that "' c '" uses IEEE floating point ops.']);
	end;
end;
return;
%_______________________________________________________________________

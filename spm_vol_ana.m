function V = spm_vol_ana(fname, n)
% Get header information etc for Analyze 7.5 image
%
% FORMAT V = spm_vol_ana(fname, n)
% fname - name of Analyze file
% n     - frame number (defaults to 1)
% V     -  a vector of structures containing image volume information.
%          See spm_vol.m for more information
%
% Note that spm_flip_analyze_images.m should be configured for your
% data.  Images that were flipped at spatial normalisation in SPM99
% should be set so that they are left-right flipped.  Images that
% were not flipped at spatial normalisation should be set to be
% unflipped.
%
% Most fields in an Analyze .hdr are ignored, except...
%
% 	* hdr.dime.datatype is used to determine the datatype of the file.
% 	* hdr.dime.dim(1)<0 | hdr.dime.dim(1)>15 indicate that values in the
% 	  .hdr and .img should be flipped.
% 	* hdr.dime.dim(2) - hdr.dime.dim(5) are the x,y,z and t dimensions.
% 	* If the .hdr is not customised for storing orientation information
% 	  then the following fields may be read to derive the information:
% 		* .mat file if it exists
% 		* otherwise hdr.dime.dim(2:4) to get origin in centre of volume
% 		  and hdr.dime.pixdim(2:4) to get the voxel sizes.
% 	* Scalefactors are derived from hdr.dime.funused1 if it is set,
% 	  otherwise they are derived from hdr.dime.cal_max, hdr.dime.cal_min,
% 	  hdr.dime.glmax and hdr.dime.glmin.
% 	* hdr.dime.vox_offset to get the offset into the volume.
% 	* hdr.hk.sizeof_hdr > 348 is used to determine if orientation
% 	  information is stored.
%
%_______________________________________________________________________
% %W% John Ashburner (from Karl's old code) %E%

if nargin<2, n = 1; end;

% Ensure that suffix is .hdr
%-----------------------------------------------------------------------
fname = deblank(fname);
q     = length(fname);
if q>=4 & fname(q - 3) == '.'; fname = fname(1:(q - 4)); end
fname = [fname '.hdr'];

% Read the .hdr information
%-----------------------------------------------------------------------
[hdr,otherendian] = read_hdr(fname);

if hdr.dime.dim(5)<n,
	error(['Not enough volumes in "' fname '" (' num2str(n) '>' num2str(hdr.dime.dim(5)) ').']);
end;

% Datatype - set image to be byte-swapped if necessary
%-----------------------------------------------------------------------
dt  = hdr.dime.datatype;
if otherendian & dt~=2 & dt~=2+128, dt = dt*256; end;
dim     = [hdr.dime.dim(2:4)' dt];
%-----------------------------------------------------------------------

% Orientation information...
%-----------------------------------------------------------------------
if hdr.hist.orient, warning(['Ignoring hdr.hist.orient field of "' fname '".']); end;

if isfield(hdr,'spmf'),
	mat = hdr.spmf(n).mat;
else,
	if exist([spm_str_manip(fname,'sd') '.mat'])==2,
		t = load([spm_str_manip(fname,'sd') '.mat']);
		mat = t.M;
	else,
		origin = (hdr.dime.dim(2:4)+1)/2;
		vox    = hdr.dime.pixdim(2:4);
		if all(vox == 0), vox = [1 1 1]; end;
		off    = -vox.*origin;
		mat    = [vox(1) 0 0 off(1) ; 0 vox(2) 0 off(2) ; 0 0 vox(3) off(3) ; 0 0 0 1];
	end;
	if spm_flip_analyze_images, mat = diag([-1 1 1 1])*mat; end;
end;
%-----------------------------------------------------------------------

% Scaling etc... and offset into file...
%-----------------------------------------------------------------------
if hdr.dime.funused1,
	scal  = hdr.dime.funused1;
	dcoff = 0;
else
	if hdr.dime.glmax-hdr.dime.glmin & hdr.dime.cal_max-hdr.dime.cal_min,
		scal  = (hdr.dime.cal_max-hdr.dime.cal_min)/(hdr.dime.glmax-hdr.dime.glmin);
		dcoff = hdr.dime.cal_min - scal*hdr.dime.glmin;
		if hdr.dime.funused1,
			warning(['Taking scalefactor etc from glmax/glmin & cal_max/cal_min of "' fname '".']);
		end;
	else,
		scal  = 1;
		dcoff = 0;
		warning(['Assuming a scalefactor of 1 for "' fname '".']);
	end;
end;

offset = hdr.dime.vox_offset + (n-1)*spm_type(hdr.dime.datatype,'bits')/8;

if rem(offset, spm_type(hdr.dime.datatype,'bits')/8),
	error(['Offset into file must be a multiple of ' spm_type(hdr.dime.datatype,'bits')/8 '.']);
end;

pinfo   = [scal dcoff offset]';
%-----------------------------------------------------------------------

% Make volume handle
%-----------------------------------------------------------------------
fname   = [spm_str_manip(fname,'sd') '.img'];
descrip = hdr.hist.descrip;
V       = struct(...
	'fname',	fname,...
	'dim',		dim,...
	'pinfo',	pinfo,...
	'mat',		mat,...
	'descrip',	descrip,...
	'n',		n,...
	'hdr',		hdr);
%-----------------------------------------------------------------------
return;
%_______________________________________________________________________
%_______________________________________________________________________
function [hdr,otherendian] = read_hdr(fname)
fid   = fopen(fname,'r','native');
otherendian = 0;
if (fid > 0)
	dime = read_dime(fid);
	if dime.dim(1)<0 | dime.dim(1)>15, % Appears to be other-endian
		% Re-open other-endian
		fclose(fid);
		if spm_platform('bigend'), fid = fopen(fname,'r','ieee-le');
		else,                      fid = fopen(fname,'r','ieee-be'); end;
		otherendian = 1;
		dime = read_dime(fid);
	end;
	hk       = read_hk(fid);
	hist     = read_hist(fid);
	hdr.hk   = hk;
	hdr.dime = dime;
	hdr.hist = hist;

	% SPM specific bit
	if hdr.hk.sizeof_hdr > 348,
		spmf = read_spmf(fid,dime.dim(5));
		if ~isempty(spmf),
			hdr.spmf = spmf;
		end;
	end;

	fclose(fid);
else,
	error(['Problem opening header file (' fopen(fid) ').']);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function hk = read_hk(fid)
% read (struct) header_key
%-----------------------------------------------------------------------
fseek(fid,0,'bof');
hk.sizeof_hdr 		= fread(fid,1,'int32');
hk.data_type  		= mysetstr(fread(fid,10,'uchar'))';
hk.db_name    		= mysetstr(fread(fid,18,'uchar'))';
hk.extents    		= fread(fid,1,'int32');
hk.session_error	= fread(fid,1,'int16');
hk.regular			= mysetstr(fread(fid,1,'uchar'))';
hk.hkey_un0			= mysetstr(fread(fid,1,'uchar'))';
if isempty(hk.hkey_un0), error(['Problem reading "hk" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function dime = read_dime(fid)
% read (struct) image_dimension
%-----------------------------------------------------------------------
fseek(fid,40,'bof');
dime.dim		= fread(fid,8,'int16');
dime.vox_units	= mysetstr(fread(fid,4,'uchar'))';
dime.cal_units	= mysetstr(fread(fid,8,'uchar'))';
dime.unused1	= fread(fid,1,'int16');
dime.datatype	= fread(fid,1,'int16');
dime.bitpix		= fread(fid,1,'int16');
dime.dim_un0	= fread(fid,1,'int16');
dime.pixdim		= fread(fid,8,'float');
dime.vox_offset	= fread(fid,1,'float');
dime.funused1	= fread(fid,1,'float');
dime.funused2	= fread(fid,1,'float');
dime.funused3	= fread(fid,1,'float');
dime.cal_max	= fread(fid,1,'float');
dime.cal_min	= fread(fid,1,'float');
dime.compressed	= fread(fid,1,'int32');
dime.verified	= fread(fid,1,'int32');
dime.glmax		= fread(fid,1,'int32');
dime.glmin		= fread(fid,1,'int32');
if isempty(dime.glmin), error(['Problem reading "dime" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function hist = read_hist(fid)
% read (struct) data_history
%-----------------------------------------------------------------------
fseek(fid,148,'bof');
hist.descrip	= mysetstr(fread(fid,80,'uchar'))';
hist.aux_file	= mysetstr(fread(fid,24,'uchar'))';
hist.orient		= fread(fid,1,'uchar');
hist.origin		= fread(fid,5,'int16');
hist.generated	= mysetstr(fread(fid,10,'uchar'))';
hist.scannum	= mysetstr(fread(fid,10,'uchar'))';
hist.patient_id	= mysetstr(fread(fid,10,'uchar'))';
hist.exp_date	= mysetstr(fread(fid,10,'uchar'))';
hist.exp_time	= mysetstr(fread(fid,10,'uchar'))';
hist.hist_un0	= mysetstr(fread(fid,3,'uchar'))';
hist.views		= fread(fid,1,'int32');
hist.vols_added	= fread(fid,1,'int32');
hist.start_field= fread(fid,1,'int32');
hist.field_skip	= fread(fid,1,'int32');
hist.omax		= fread(fid,1,'int32');
hist.omin		= fread(fid,1,'int32');
hist.smax		= fread(fid,1,'int32');
hist.smin		= fread(fid,1,'int32');
if isempty(hist.smin), error(['Problem reading "hist" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function spmf = read_spmf(fid,n)
% Read SPM specific fields
fseek(fid,348,'bof');
mgc = fread(fid,1,'int32');    % Magic number
if mgc ~= 20020417, spm = []; return; end;

for j=1:n,
	spmf(j).mat    = fread(fid,16,'double'); % Orientation information
	spmf(j).unused = fread(fid,384,'uchar'); % Extra unused stuff
	if length(spmf(j).unused)<384,
		error(['Problem reading "spmf" of header file (' fopen(fid) ').']);
	end;
 	spmf(j).mat = reshape(spmf(j).mat,[4 4]);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function out = mysetstr(in)
tmp = find(in == 0);
tmp = min([min(tmp) length(in)]);
out = setstr(in(1:tmp));
return;
%_______________________________________________________________________
%_______________________________________________________________________

function [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = hread(P)
% reads a header
% FORMAT [DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P);
%
% P       - filename 	     (e.g spm or spm.img)
% DIM     - image size       [i j k [l]] (voxels)
% VOX     - voxel size       [x y z [t]] (mm [secs])
% SCALE   - scale factor
% TYPE    - datatype (integer - see spm_type)
% OFFSET  - offset (bytes)
% ORIGIN  - origin [i j k]
% DESCRIP - description string
%___________________________________________________________________________
%
% spm_hread reads variables into working memory from a SPM/ANALYZE
% compatible header file.  If the header does not exist global defaults
% are used.  The 'originator' field of the ANALYZE format has been
% changed to ORIGIN in the SPM version of the header.  funused1 of the
% ANALYZE format is used for SCALE
%
% see also dbh.h (ANALYZE) spm_hwrite.m and spm_type.m
%
%__________________________________________________________________________
% %W% %E%

return;
%_______________________________________________________________________

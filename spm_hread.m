function [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = spm_hread(P)
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

% ensure correct suffix {.hdr}
%---------------------------------------------------------------------------
P     = P(P ~= ' ');
q     = length(P);
if P(q - 3) == '.'; P = P(1:(q - 4)); end
P     = [P '.hdr'];

% open header file
%---------------------------------------------------------------------------
P     = P(P ~= ' ');
fid   = fopen(P,'r','native');

if (fid > 0)

% read (struct) header_key
%---------------------------------------------------------------------------
fseek(fid,0,'bof');

otherendian = 0;
sizeof_hdr 	= fread(fid,1,'int32');
if sizeof_hdr==1543569408, % Appears to be other-endian
	% Re-open other-endian
	fclose(fid);
	if spm_bigend,
		fid = fopen(P,'r','ieee-le');
	else,
		fid = fopen(P,'r','ieee-be');
	end;
	fseek(fid,0,'bof');
	sizeof_hdr = fread(fid,1,'int32');
	otherendian = 1;
end;

data_type  	= setstr(fread(fid,10,'char'))';
db_name    	= setstr(fread(fid,18,'char'))';
extents    	= fread(fid,1,'int32');
session_error   = fread(fid,1,'int16');
regular    	= setstr(fread(fid,1,'char'))';
hkey_un0    	= setstr(fread(fid,1,'char'))';



% read (struct) image_dimension
%---------------------------------------------------------------------------
fseek(fid,40,'bof');

dim    		= fread(fid,8,'int16');
vox_units    	= setstr(fread(fid,4,'char'))';
cal_units    	= setstr(fread(fid,8,'char'))';
unused1		= fread(fid,1,'int16');
datatype	= fread(fid,1,'int16');
bitpix		= fread(fid,1,'int16');
dim_un0		= fread(fid,1,'int16');
pixdim		= fread(fid,8,'float');
vox_offset	= fread(fid,1,'float');
funused1	= fread(fid,1,'float');
funused2	= fread(fid,1,'float');
funused3	= fread(fid,1,'float');
cal_max		= fread(fid,1,'float');
cal_min		= fread(fid,1,'float');
compressed	= fread(fid,1,'int32');
verified	= fread(fid,1,'int32');
glmax		= fread(fid,1,'int32');
glmin		= fread(fid,1,'int32');

% read (struct) data_history
%---------------------------------------------------------------------------
fseek(fid,148,'bof');

descrip		= setstr(fread(fid,80,'char'))';
aux_file	= setstr(fread(fid,24,'char'))';
orient		= fread(fid,1,'char');
origin		= fread(fid,5,'int16');
generated	= setstr(fread(fid,10,'char'))';
scannum		= setstr(fread(fid,10,'char'))';
patient_id	= setstr(fread(fid,10,'char'))';
exp_date	= setstr(fread(fid,10,'char'))';
exp_time	= setstr(fread(fid,10,'char'))';
hist_un0	= setstr(fread(fid,3,'char'))';
views		= fread(fid,1,'int32');
vols_added	= fread(fid,1,'int32');
start_field	= fread(fid,1,'int32');
field_skip	= fread(fid,1,'int32');
omax		= fread(fid,1,'int32');
omin		= fread(fid,1,'int32');
smax		= fread(fid,1,'int32');
smin		= fread(fid,1,'int32');

fclose(fid);

if isempty(smin)
	error(['There is a problem with the header file ' P '.']);
end

% convert to SPM global variables
%---------------------------------------------------------------------------
DIM    	  	= dim(2:4)';
VOX       	= pixdim(2:4)';
SCALE     	= funused1;
SCALE    	= ~SCALE + SCALE;
TYPE     	= datatype;
if otherendian == 1 & datatype ~= 2,
	TYPE = TYPE*256;
end;
OFFSET	  	= vox_offset;
ORIGIN    	= origin(1:3)';
DESCRIP   	= descrip(1:max(find(descrip)));

else
	global DIM VOX SCALE TYPE OFFSET ORIGIN
	DESCRIP = ['defaults'];
end
return;
%_______________________________________________________________________

function VO = spm_write_plane(V,A,p)
% Write a transverse plane of image data.
% FORMAT spm_write_plane(V,A,p)
% V   - data structure containing image information.
%       - see spm_vol for a description.
% A   - the two dimensional image to write.
% p   - the plane number (beginning from 1).
%
% VO  - (possibly) modified data structure containing image information.
%       It is possible that future versions of spm_write_plane may
%       modify scalefactors (for example).
%
%_______________________________________________________________________
% %W% John Ashburner %E%

if any(V.dim(1:2) ~= size(A))
	error('Incompatible image dimensions');
end;

% Write Analyze image by default
VO = write_analyze_plane(V,A,p);
if isempty(VO),
	spm_progress_bar('Clear');
	write_error_message(V.fname);
	error(['Error writing ' V.fname '.']);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function VO = write_analyze_plane(V,A,p)
VO = V;
% Check datatype is OK
dt = V.dim(4);

% Convert to native datatype
dt = spm_type(spm_type(dt));
s  = find(dt == spm_type);
if isempty(s)
	sts = -1;
	disp(['Unrecognised data type (' num2str(V.dim(4)) ')']);
	return;
end;
datasize = spm_type(dt,'bits')/8;

% Open the image file
fname = [spm_str_manip(V.fname,'svd') '.img'];
fp    = fopen(fname,'r+');
if fp == -1,
	fp = fopen(fname,'w');
	if fp == -1,
		VO = [];
		disp('Can''t open image file');
		return;
	end;
end;

% Rescale to fit appropriate range
if spm_type(dt,'intt'),
	maxval = max(spm_type(dt,'maxval')*V.pinfo(1,:) + V.pinfo(2,:));
	if any(dt == [128+2 128+4 128+8]),
		 % Convert to a form that Analyze will support
		dt = dt - 128; 
	end;
	mxv = spm_type(dt,'maxval');
	mnv = spm_type(dt,'minval');
	scale  = maxval/mxv;
	A      = round(A/scale);
	A(find(A > mxv)) = mxv;
	A(find(A < mnv)) = mnv;
else,
	scale = max(V.pinfo(1,:));
	A     = A/scale;
end;

% Seek to the appropriate offset
off   = (p-1)*datasize*prod(V.dim(1:2)); % + V.pinfo(3,1);
if fseek(fp,off,'bof')==-1,
	% Need this because fseek in Matlab does not
	% seek past the EOF
	fseek(fp,0,'eof');
	curr_off = ftell(fp);
	blanks   = zeros(off-curr_off,1);
	if fwrite(fp,blanks,'uchar') ~= prod(size(blanks)),
		fclose(fp);
		VO = [];
		return;
	end;
	if fseek(fp,off,'bof') == -1,
		fclose(fp);
		VO = [];
		return;
	end;
end;

if fwrite(fp,A,spm_type(dt)) ~= prod(size(A)),
	fclose(fp);
	VO = [];
	return;
end;
fclose(fp);

return;
%_______________________________________________________________________

%_______________________________________________________________________
function write_error_message(q)

str = {...
	'Error writing:',...
	' ',...
	['        ',spm_str_manip(q,'k40d')],...
	' ',...
	'Check file permissions / disk space / disk quota.'};
spm('alert*',str,mfilename,sqrt(-1));

return;
%_______________________________________________________________________

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
	error(['Error writing ' V.fname '. Check your disk space.']);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function VO = write_analyze_plane(V,A,p)
VO = V;
% Check datatype is OK
dt = V.dim(4);

% Convert to native datatype
if dt>256,
	dt = dt/256;
end;
s = find(dt == [2 4 8 16 64 128+2 128+4 128+8]);
if isempty(s)
	VO = [];
	disp(['Unrecognised data type (' num2str(V.dim(4)) ')']);
	return;
end;
datasize = [1 2 4 4 8 1 2 4];
datasize = datasize(s);

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
if any(dt == [2 4 8 128+2 128+4 128+8]),
	msk = find(~finite(A));
	A(msk) = 0.0;
	if dt == 2,
		maxval = max(255*V.pinfo(1,:) + V.pinfo(2,:));
		scale  = maxval/255;
		A      = round(A/scale);
		A(find(A <   0)) =   0;
		A(find(A > 255)) = 255;
	elseif dt == 4,
		maxval = max(32767*V.pinfo(1,:) + V.pinfo(2,:));
		scale  = maxval/32767;
		A      = round(A/scale);
		A(find(A >  32767)) =  32767;
		A(find(A < -32768)) = -32768;
	elseif (dt == 8)
		maxval = max((2^31-1)*V.pinfo(1,:) + V.pinfo(2,:));
		scale  = maxval/(2^31-1);
		A      = round(A/scale);
		A(find(A >  2^31-1)) =  2^31-1;
		A(find(A < -2^31  )) = -2^31;
	elseif dt == 128+2,
		maxval = max(127*V.pinfo(1,:) + V.pinfo(2,:));
		scale  = maxval/255;
		A      = round(A/scale);
		A(find(A <   0)) =   0;
		A(find(A > 255)) = 255;
		dt     = dt - 128;
	elseif dt == 128+4,
		maxval = max(65535*V.pinfo(1,:) + V.pinfo(2,:));
		scale  = maxval/32767;
		A      = round(A/scale);
		A(find(A >  32767)) =  32767;
		A(find(A < -32768)) = -32768;
		dt     = dt - 128;
	elseif (dt == 128+8)
		maxval = max((2^32-1)*V.pinfo(1,:) + V.pinfo(2,:));
		scale  = maxval/(2^31-1);
		A      = round(A/scale);
		A(find(A >  2^31-1)) =  2^31-1;
		A(find(A < -2^31  )) = -2^31;
		dt     = dt - 128;
	end;
end;

% Seek to the appropriate offset
off   = (p-1)*datasize*prod(V.dim(1:2)); % + V.pinfo(3,1);
if fseek(fp,off,'bof')==-1,
	% Need this because fseek in Matlab does not
	% seek past the EOF
	fseek(fp,0,'eof');
	curr_off = ftell(fp);
	blanks = zeros(off-curr_off,1);
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
	'Error opening:',...
	['        ',spm_str_manip(q,'k40d')],...
	' ',...
	'Check file permissions / disk space / disk quota.'};

msgbox(str,sprintf('%s%s: %s...',spm('ver'),spm('GetUser',' (%s)'),mfilename),...
	'error')

return;
%_______________________________________________________________________

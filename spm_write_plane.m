function spm_write_plane(V,A,p)
% Write a transverse plane of image data.
% FORMAT spm_write_plane(V,A,p)
% V   - data structure containing image information.
%       - see spm_vol for a description.
% A   - the two dimensional image to write.
% p   - the plane number (beginning from 1).
%
%_______________________________________________________________________
% %W% John Ashburner %E%

if any(V.dim(1:2) ~= size(A))
	error('Incompatible image dimensions');
end;

% Write Analyze image by default
if write_analyze_plane(V,A,p) == -1,
	spm_progress_bar('Clear');
	write_error_message(V.fname);
	error(['Error writing ' V.fname '. Check your disk space.']);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function sts = write_analyze_plane(V,A,p)
sts   = 0;

% Check datatype is OK
s = find(V.dim(4) == [2 4 8 16 64]);
if isempty(s)
	sts = -1;
	disp(['Unrecognised data type (' num2str(V.dim(4)) ')']);
	return;
end;
datasize = [1 2 4 4 8];
datasize = datasize(s);

% Open the image file
fname = [spm_str_manip(V.fname,'svd') '.img'];
fp    = fopen(fname,'r+');
if fp == -1,
	fp = fopen(fname,'w');
	if fp == -1,
		sts = -1;
		disp('Can''t open image file');
		return;
	end;
end;

% Rescale for writing
A = A/V.pinfo(1,1);

% Truncate and round if necessary
if any(V.dim(4) == [2 4 8]),
	A = round(A);
	if V.dim(4) == 2,
		A(find(A <   0)) =   0;
		A(find(A > 255)) = 255;
	elseif (V.dim(4) == 4)
		A(find(A >  32767)) =  32767;
		A(find(A < -32768)) = -32768;
	elseif (V.dim(4) == 8)
		A(find(A >  2^31-1)) =  2^31-1;
		A(find(A < -2^31  )) = -2^31;
	end
end;

% Write the plane at the appropriate offset
off   = V.pinfo(3,1) + (p-1)*datasize*prod(V.dim(1:2));
fseek(fp,off,'bof');
if fwrite(fp,A,spm_type(V.dim(4))) ~= prod(size(A)),
	fclose(fp);
	sts = -1;
	return;
end;
fclose(fp);

return;
%_______________________________________________________________________

%_______________________________________________________________________
function write_error_message(q)
f=spm_figure('findwin','Graphics'); 
if ~isempty(f), 
	figure(f); 
	spm_figure('Clear','Graphics'); 
	spm_figure('Clear','Interactive'); 
	ax=axes('Visible','off','Parent',f); 
	text(0,0.60,'Error opening:', 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.55,spm_str_manip(q,'k40d'), 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.40,'  Please check that you have write permission.', 'FontSize', 16, 'Interpreter', 'none'); 
end
return;
%_______________________________________________________________________

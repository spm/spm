function spm_create_image(V)
% Create an image file.
% FORMAT spm_create_image(V)
% V   - data structure containing image information.
%       - see spm_vol for a description.
%
%_______________________________________________________________________
% %W% John Ashburner %E%

if create_analyze_image(V) == -1,
	spm_progress_bar('Clear');
	open_error_message(V.fname);
	error(['Error opening ' V.fname '. Check that you have write permission.']);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function sts = create_analyze_image(V)
sts  = 0;

% Extract voxel sizes and origin from the 4x4 matrix
vx   = sqrt(sum(V.mat(1:3,1:3).^2));
orgn = V.mat\[0 0 0 1]';
orgn = orgn(1:3)';

% Add a description when available
if isfield(V,'descrip'),
	descrip = V.descrip;
else
	descrip = 'SPM compatible';
end;
% Check datatype is OK
dt = V.dim(4);

% Convert to native datatype
if dt>256,
	dt = dt/256;
end;
s = find(dt == [2 4 8 16 64 128+2 128+4 128+8]);
if isempty(s)
	sts = -1;
	disp(['Unrecognised data type (' num2str(V.dim(4)) ')']);
	return;
end;

% Compute an appropriate scalefactor
if dt == 2,
	maxval = max(255*V.pinfo(1,:) + V.pinfo(2,:));
	scale  = maxval/255;
elseif dt == 4,
	maxval = max(32767*V.pinfo(1,:) + V.pinfo(2,:));
	scale  = maxval/32767;
elseif (dt == 8)
	maxval = max((2^31-1)*V.pinfo(1,:) + V.pinfo(2,:));
	scale  = maxval/(2^31-1);
elseif dt == 128+2,
	maxval = max(127*V.pinfo(1,:) + V.pinfo(2,:));
	scale  = maxval/255;
	dt     = dt - 128;
elseif dt == 128+4,
	maxval = max(65535*V.pinfo(1,:) + V.pinfo(2,:));
	scale  = maxval/32767;
	dt     = dt - 128;
elseif (dt == 128+8)
	maxval = max((2^32-1)*V.pinfo(1,:) + V.pinfo(2,:));
	scale  = maxval/(2^31-1);
	dt     = dt - 128;
else
	scale = 1.0;
end;


% Write the header
s    = spm_hwrite(deblank(V.fname), [V.dim(1:3) 1],...
	vx, scale, dt, 0, orgn, descrip);
if s~= 348,
	sts = -1;
end;

% Write the matrix
spm_get_space(deblank(V.fname),V.mat);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function open_error_message(q)
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
return
%_______________________________________________________________________

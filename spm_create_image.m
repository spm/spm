function Vo=spm_create_image(Vi)
% Create an image file.
% FORMAT Vo = spm_create_image(Vi)
% Vi   - data structure containing image information.
%      - see spm_vol for a description.
% Vo   - data structure after modification for writing.
%_______________________________________________________________________
% %W% John Ashburner %E%

[sts,Vo] =  create_analyze_image(Vi);
if sts == -1,
	spm_progress_bar('Clear');
	open_error_message(Vi.fname);
	error(['Error opening ' Vi.fname '. Check that you have write permission.']);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [sts,V] = create_analyze_image(V)
sts  = 0;

% Extract voxel sizes and origin from the 4x4 matrix
vx   = sqrt(sum(V.mat(1:3,1:3).^2));
orgn = V.mat\[0 0 0 1]';
orgn = orgn(1:3)';

% Add a description when available
if ~isfield(V,'descrip'),
	V.descrip = 'SPM compatible';
end;
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

% Compute an appropriate scalefactor
if spm_type(dt,'intt'),
	maxval = max(spm_type(dt,'maxval')*V.pinfo(1,:) + V.pinfo(2,:));
	if any(dt == [128+2 128+4 128+8]),
		 % Convert to a form that Analyze will support
		dt = dt - 128; 
	end;
	scale  = maxval/spm_type(dt,'maxval');
else
	scale = max(V.pinfo(1,:));
end;

V.pinfo = [scale 0 0]';
V.dim(4)= dt;

% Write the header
s    = spm_hwrite(deblank(V.fname), [V.dim(1:3) 1],...
	vx, scale, dt, 0, orgn, V.descrip);
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

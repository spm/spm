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

% Write the header
s    = spm_hwrite(deblank(V.fname), [V.dim(1:3) 1],...
	vx, V.pinfo(1,1), V.dim(4), V.pinfo(3,1), orgn, descrip);
if s~= 348,
	sts = -1;
end;

% Write the matrix
spm_get_space(deblank(V.fname),V.mat);
return;

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

function V = spm_vol(P)
% Get header information etc for images.
% FORMAT V = spm_vol(P)
% P - a matrix of filenames.
% V - a vector of structures containing image volume information.
% The elements of the structures are:
%       V.fname - the filename of the image.
%       V.dim   - the x, y and z dimensions of the volume, and the
%                 datatype of the image.
%       V.mat   - a 4x4 affine transformation matrix mapping from
%                 voxel coordinates to real world coordinates.
%       V.pinfo - plane info for each plane of the volume.
%              V.pinfo(1,:) - scale for each plane
%              V.pinfo(2,:) - offset for each plane
%                 The true voxel intensities of the jth image are given
%                 by: val*V.pinfo(1,j) + V.pinfo(2,j)
%              V.pinfo(3,:) - offset into image (in bytes).
%                 If the size of pinfo is 3x1, then the volume is assumed
%                 to be contiguous and each plane has the same scalefactor
%                 and offset.
%____________________________________________________________________________
%
% The fields listed above are essential for the mex routines, but other
% fields can also be incorporated into the structure.
% Unlike with spm_map_vol, the images are not memory mapped at this step,
% but are mapped when the mex routines using the volume information are
% called.
% This is a replacement for the spm_map_vol and spm_unmap_vol stuff
% that I expect will become redundant.
%_______________________________________________________________________
% %W% John Ashburner %E%

V = subfunc2(P);
return;

function V = subfunc2(P)
if iscell(P),
	V = cell(size(P));
	for j=1:prod(size(P)),
		if iscell(P{j}),
			V{j} = subfunc2(P{j});
		else,
			V{j} = subfunc1(P{j});
		end;
	end;
else
	V = subfunc1(P);
end;
return;

function V = subfunc1(P)
if size(P,1)==0,
	V=[];
else,
	V(size(P,1),1) = struct('fname','', 'dim', [0 0 0 0], 'mat',eye(4), 'pinfo', [1 0 0]');
end;
for i=1:size(P,1),
	v = subfunc(P(i,:));
	if isempty(v),
		hread_error_message(P(i,:));
		error(['Can''t get volume information for ''' P(i,:) '''']);
	end;
	f = fieldnames(v);
	for j=1:size(f,1),
		eval(['V(i).' f{j} ' = v.' f{j} ';']);
	end;
end;
return;

function V = subfunc(p)
p = deblank(p);

if exist([spm_str_manip(p,'sd') '.hdr']) == 2,
	[dim vox scale dtype offset origin descrip] = spm_hread(p);
	mat = spm_get_space([spm_str_manip(p,'sd') '.img']);
	V   = struct('fname',[spm_str_manip(p,'sd') '.img'],...
	           'dim',[dim dtype],'mat',mat,'pinfo',[scale 0 offset]',...
	           'descrip',descrip);

else, % Try other formats

	% Try MINC format
	V=spm_vol_minc(p);
	if ~isempty(V), return; end;

	% Try Ecat 7
	V=spm_vol_ecat7(p);
	if ~isempty(V), return; end;
end;
return;


%_______________________________________________________________________
function hread_error_message(q)
f=spm_figure('findwin','Graphics'); 
if ~isempty(f), 
	figure(f); 
	spm_figure('Clear','Graphics'); 
	spm_figure('Clear','Interactive'); 
	ax=axes('Visible','off','Parent',f); 
	text(0,0.60,'Error reading information on:', 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.55,spm_str_manip(q,'k40d'), 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.40,'  Please check that it is in the correct format.', 'FontSize', 16, 'Interpreter', 'none'); 
end
return;
%_______________________________________________________________________

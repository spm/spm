function spm_image(op,varargin)
% image and header display
% FORMAT spm_image
%_______________________________________________________________________
%
% spm_image is an interactive facility that allows orthogonal sections
% from an image volume to be displayed.  Clicking the cursor on either
% of the three images moves the point around which the orthogonal
% sections are viewed.  The co-ordinates of the cursor are shown both
% in voxel co-ordinates and millimeters within some fixed framework.
% The intensity at that point in the image (sampled using tri-linear
% interpolation) is also given. The position of the crosshairs can also
% be moved by specifying the co-ordinates in millimeters to which they
% should be moved.
%
% The images can be re-oriented by entering appropriate translations,
% rotations and zooms into the panel on the left.  The transformations
% can then be saved by hitting the ``Reorient images...'' button.  The
% transformations that were applied to the image are saved to the
% ``.mat'' files of the selected images.  The transformations are
% considered to be relative to any existing transformations that may be
% stored in the ``.mat'' files.  Note that the order that the
% transformations are applied in is the same as in ``spm_matrix.m''.
%
% The right panel shows miscellaneous information about the image.
% This includes:
%   Dimensions - the x, y and z dimensions of the image.
%   Datatype   - the computer representation of each voxel.
%   Intensity  - scalefactors and possibly a DC offset.
%   Miscellaneous other information about the image.
%   Vox size   - the distance (in mm) between the centres of
%                neighbouring voxels.
%   Origin     - the voxel at the origin of the co-ordinate system
%   DIr Cos    - Direction cosines.  This is a widely used
%                representation of the orientation of an image.
%
% There are also a few options for different resampling modes, zooms
% etc.  You can also flip between voxel space (as would be displayed
% by Analyze) or world space (the orientation that SPM considers the
% image to be in).  If you are re-orienting the images, make sure that
% world space is specified.  Blobs (from activation studies) can be
% superimposed on the images and the intensity windowing can also be
% changed.
%
%_______________________________________________________________________
% %W% John Ashburner %E%

global st

if nargin==0,
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Display',0);
	SPMid = spm('FnBanner',mfilename,'%I%');
	spm_help('!ContextHelp',[mfilename,'.m'])

	% get the image's filename {P}
	%-----------------------------------------------------------------------
	P      = spm_get(1,'.img','please select image',[]);
	spm_image('init',P);
	return;
end;

if ~strcmp(op,'init') & ~strcmp(op,'reset') & isempty(st.vols{1}), my_reset; warning('Lost all the image information'); return; end;

if strcmp(op,'repos'),
	% The widgets for translation rotation or zooms have been modified.
	%-----------------------------------------------------------------------
	fg = spm_figure('Findwin','Graphics');
	set(fg,'Pointer','watch');
	i = varargin{1};
	st.B(i) = eval(get(gco,'String'),num2str(st.B(i)));
	set(gco,'String',st.B(i));
	st.vols{1}.premul = spm_matrix(st.B);
	% spm_orthviews('MaxBB');
	spm_image('zoom_in');
	spm_image('update_info');
	set(fg,'Pointer','arrow');
	return;
end;

if strcmp(op,'shopos'),
	% The position of the crosshairs has been moved.
	%-----------------------------------------------------------------------
	if isfield(st,'mp'),
		fg = spm_figure('Findwin','Graphics');
		if any(findobj(fg) == st.mp),
		set(st.mp,'String',sprintf('%.1f %.1f %.1f',spm_orthviews('pos')));
		pos = spm_orthviews('pos',1);
		set(st.vp,'String',sprintf('%.1f %.1f %.1f',pos));
		set(st.in,'String',sprintf('%g',spm_sample_vol(st.vols{1},pos(1),pos(2),pos(3),st.hld)));
		else,
			st.Callback = ';';
			rmfield(st,{'mp','vp','in'});
		end;
	else,
		st.Callback = ';';
	end;
	return;
end;

if strcmp(op,'setposmm'),
	% Move the crosshairs to the specified position
	%-----------------------------------------------------------------------
	if isfield(st,'mp'),
		fg = spm_figure('Findwin','Graphics');
		if any(findobj(fg) == st.mp),
			pos = sscanf(get(st.mp,'String'), '%g %g %g');
			if length(pos)~=3,
				pos = spm_orthviews('pos');
			end;
			spm_orthviews('Reposition',pos);
		end;
	end;
	return;
end;

if strcmp(op,'setposvx'),
	% Move the crosshairs to the specified position
	%-----------------------------------------------------------------------
	if isfield(st,'mp'),
		fg = spm_figure('Findwin','Graphics');
		if any(findobj(fg) == st.vp),
			pos = sscanf(get(st.vp,'String'), '%g %g %g');
			if length(pos)~=3,
				pos = spm_orthviews('pos',1);
			end;
			tmp = st.vols{1}.premul*st.vols{1}.mat;
			pos = tmp(1:3,:)*[pos ; 1];
			spm_orthviews('Reposition',pos);
		end;
	end;
	return;
end;


if strcmp(op,'addblobs'),
	% Add blobs to the image - in full colour
	spm_figure('Clear','Interactive');
	nblobs = spm_input('Number of sets of blobs',1,'1|2|3|4|5|6',[1 2 3 4 5 6],1);
	for i=1:nblobs,
		[SPM,VOL,DES] = spm_getSPM;
		c = spm_input('Colour','+1','m','Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
		colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
		spm_orthviews('addcolouredblobs',1,SPM.XYZ,SPM.Z,VOL.M,colours(c,:));
		set(st.blobber,'String','Remove Blobs','Callback','spm_image(''rmblobs'');');
	end;
	spm_orthviews('Redraw');
end;

if strcmp(op,'rmblobs'),
	% Remove all blobs from the images
	spm_orthviews('rmblobs',1);
	set(st.blobber,'String','Add Blobs','Callback','spm_image(''addblobs'');');
end;

if strcmp(op,'window'),
	op = get(st.win,'Value');
	if op == 1,
		spm_orthviews('window',1);
	else,
		spm_orthviews('window',1,spm_input('Range','+1','e','',2));
	end;
end;


if strcmp(op,'reorient'),
	% Time to modify the ``.mat'' files for the images.
	% I hope that giving people this facility is the right thing to do....
	%-----------------------------------------------------------------------
	mat = spm_matrix(st.B);
	if det(mat)<=0,
		h=msgbox('This will flip the images','Warning!','warn');
		uiwait(h);
	end;
	P = spm_get(Inf, 'img','Images to reorient');
	Mats = zeros(4,4,size(P,1));
	spm_progress_bar('Init',size(P,1),'Reading current orientations',...
		'Images Complete');
	for i=1:size(P,1),
		Mats(:,:,i) = spm_get_space(P(i,:));
		spm_progress_bar('Set',i);
	end;
	spm_progress_bar('Init',size(P,1),'Reorienting images',...
		'Images Complete');
	for i=1:size(P,1),
		spm_get_space(P(i,:),mat*Mats(:,:,i));
		spm_progress_bar('Set',i);
	end;
	spm_progress_bar('Clear');
	tmp = spm_get_space(st.vols{1}.fname);
	if sum((tmp(:)-st.vols{1}.mat(:)).^2) > 1e-8,
		spm_image('init',st.vols{1}.fname);
	end;
	return;
end;

if strcmp(op,'update_info'),
	% Modify the positional information in the right hand panel.
	%-----------------------------------------------------------------------
	mat = st.vols{1}.premul*st.vols{1}.mat;
	Z = spm_imatrix(mat);
	Z = Z(7:9);

	set(st.posinf.z,'String', sprintf('%.3g x %.3g x %.3g', Z));

	O = mat\[0 0 0 1]'; O=O(1:3)';
	set(st.posinf.o, 'String', sprintf('%.3g %.3g %.3g', O));

	R = spm_imatrix(mat);
	R = spm_matrix([0 0 0 R(4:6)]);
	R = R(1:3,1:3);

	tmp2 = sprintf('%+5.3f %+5.3f %+5.3f', R(1,1:3)); tmp2(find(tmp2=='+')) = ' ';
	set(st.posinf.m1, 'String', tmp2);
	tmp2 = sprintf('%+5.3f %+5.3f %+5.3f', R(2,1:3)); tmp2(find(tmp2=='+')) = ' ';
	set(st.posinf.m2, 'String', tmp2);
	tmp2 = sprintf('%+5.3f %+5.3f %+5.3f', R(3,1:3)); tmp2(find(tmp2=='+')) = ' ';
	set(st.posinf.m3, 'String', tmp2);

	tmp = [[R zeros(3,1)] ; 0 0 0 1]*diag([Z 1])*spm_matrix(-O) - mat;

	if sum(tmp(:).^2)>1e-8,
		set(st.posinf.w, 'String', 'Warning: shears involved');
	else,
		set(st.posinf.w, 'String', '');
	end;

	return;
end;

if strcmp(op,'reset'),
	my_reset;
end;

if strcmp(op,'zoom_in'),
	op = get(st.zoomer,'Value');
	if op==1,
		spm_orthviews('resolution',1);
		spm_orthviews('MaxBB');
	else,
		vx = sqrt(sum(st.Space(1:3,1:3).^2));
		vx = vx.^(-1);
		pos = spm_orthviews('pos');
		pos = st.Space\[pos ; 1];
		pos = pos(1:3)';
		if     op == 2, st.bb = [pos-80*vx ; pos+80*vx] ; spm_orthviews('resolution',1);
		elseif op == 3, st.bb = [pos-40*vx ; pos+40*vx] ; spm_orthviews('resolution',.5);
		elseif op == 4, st.bb = [pos-20*vx ; pos+20*vx] ; spm_orthviews('resolution',.25);
		elseif op == 5, st.bb = [pos-10*vx ; pos+10*vx] ; spm_orthviews('resolution',.125);
		else          , st.bb = [pos- 5*vx ; pos+ 5*vx] ; spm_orthviews('resolution',.125);
		end;
	end;
	return;
end;

if strcmp(op,'init'),
fg = spm_figure('Findwin','Graphics');
if isempty(fg)
	fg=spm_figure('Create','Graphics');
	if isempty(fg)
		error('Cant create graphics window');
	end
else,
	spm_figure('Clear','Graphics');
end;

P = deblank(varargin{1});
spm_orthviews('Reset');
spm_orthviews('Image', P, [0.0 0.45 1 0.55]);
if isempty(st.vols{1}), return; end;

spm_orthviews('MaxBB');
st.callback = 'spm_image(''shopos'');';

st.B = [0 0 0  0 0 0  1 1 1  0 0 0];

% locate Graphics window and clear it
%-----------------------------------------------------------------------
WS = spm('WinScale');

% Widgets for re-orienting images.
%-----------------------------------------------------------------------
uicontrol(fg,'Style','Frame','Position',[60 25 200 325].*WS,'DeleteFcn','spm_image(''reset'');');
uicontrol(fg,'Style','Text', 'Position',[90 220 100 016].*WS,'String','right  {mm}');
uicontrol(fg,'Style','Text', 'Position',[90 200 100 016].*WS,'String','foward  {mm}');
uicontrol(fg,'Style','Text', 'Position',[90 180 100 016].*WS,'String','up  {mm}');
uicontrol(fg,'Style','Text', 'Position',[90 160 100 016].*WS,'String','pitch  {rad}');
uicontrol(fg,'Style','Text', 'Position',[90 140 100 016].*WS,'String','roll  {rad}');
uicontrol(fg,'Style','Text', 'Position',[90 120 100 016].*WS,'String','yaw  {rad}');
uicontrol(fg,'Style','Text', 'Position',[90 100 100 016].*WS,'String','resize  {x}');
uicontrol(fg,'Style','Text', 'Position',[90  80 100 016].*WS,'String','resize  {y}');
uicontrol(fg,'Style','Text', 'Position',[90  60 100 016].*WS,'String','resize  {z}');

uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',1)','Position',[190 220 050 020].*WS,'String','0','ToolTipString','translate');
uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',2)','Position',[190 200 050 020].*WS,'String','0','ToolTipString','translate');
uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',3)','Position',[190 180 050 020].*WS,'String','0','ToolTipString','translate');
uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',4)','Position',[190 160 050 020].*WS,'String','0','ToolTipString','rotate');
uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',5)','Position',[190 140 050 020].*WS,'String','0','ToolTipString','rotate');
uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',6)','Position',[190 120 050 020].*WS,'String','0','ToolTipString','rotate');
uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',7)','Position',[190 100 050 020].*WS,'String','1','ToolTipString','zoom');
uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',8)','Position',[190  80 050 020].*WS,'String','1','ToolTipString','zoom');
uicontrol(fg,'Style','edit','Callback','spm_image(''repos'',9)','Position',[190  60 050 020].*WS,'String','1','ToolTipString','zoom');

uicontrol(fg,'Style','Pushbutton','String','Reorient images...','Callback','spm_image(''reorient'')',...
         'Position',[90 35 150 020].*WS,'ToolTipString','modify position information of selected images');

% Crosshair position
%-----------------------------------------------------------------------
uicontrol(fg,'Style','Frame','Position',[70 250 180 90].*WS);
uicontrol(fg,'Style','Text', 'Position',[75 320 170 016].*WS,'String','Crosshair Position');
uicontrol(fg,'Style','Text', 'Position',[75 300 35 016].*WS,'String','mm:');
uicontrol(fg,'Style','Text', 'Position',[75 280 35 016].*WS,'String','vx:');
uicontrol(fg,'Style','Text', 'Position',[75 260 65 016].*WS,'String','Intensity:');

st.mp = uicontrol(fg,'Style','edit', 'Position',[110 300 135 020].*WS,'String','','Callback','spm_image(''setposmm'')','ToolTipString','move crosshairs to mm coordinates');
st.vp = uicontrol(fg,'Style','edit', 'Position',[110 280 135 020].*WS,'String','','Callback','spm_image(''setposvx'')','ToolTipString','move crosshairs to voxel coordinates');
st.in = uicontrol(fg,'Style','Text', 'Position',[140 260  85 016].*WS,'String','');

% General information
%-----------------------------------------------------------------------
uicontrol(fg,'Style','Frame','Position',[305  25 280 325].*WS);
uicontrol(fg,'Style','Text','Position' ,[310 330 50 016].*WS,...
	'HorizontalAlignment','right', 'String', 'File:');
uicontrol(fg,'Style','Text','Position' ,[360 330 210 016].*WS,...
	'HorizontalAlignment','left', 'String', spm_str_manip(st.vols{1}.fname,'k25'),'FontWeight','bold');
uicontrol(fg,'Style','Text','Position' ,[310 310 100 016].*WS,...
	'HorizontalAlignment','right', 'String', 'Dimensions:');
uicontrol(fg,'Style','Text','Position' ,[410 310 160 016].*WS,...
	'HorizontalAlignment','left', 'String', sprintf('%d x %d x %d', st.vols{1}.dim(1:3)),'FontWeight','bold');
uicontrol(fg,'Style','Text','Position' ,[310 290 100 016].*WS,...
	'HorizontalAlignment','right', 'String', 'Datatype:');
uicontrol(fg,'Style','Text','Position' ,[410 290 160 016].*WS,...
	'HorizontalAlignment','left', 'String', spm_type(st.vols{1}.dim(4)),'FontWeight','bold');
uicontrol(fg,'Style','Text','Position' ,[310 270 100 016].*WS,...
	'HorizontalAlignment','right', 'String', 'Intensity:');
str = 'varied';
if size(st.vols{1}.pinfo,2) == 1,
	if st.vols{1}.pinfo(2),
		str = sprintf('Y = %g X + %g', st.vols{1}.pinfo(1:2)');
	else,
		str = sprintf('Y = %g X', st.vols{1}.pinfo(1)');
	end;
end;
uicontrol(fg,'Style','Text','Position' ,[410 270 160 016].*WS,...
	'HorizontalAlignment','left', 'String', str,'FontWeight','bold');

if isfield(st.vols{1}, 'descrip'),
	uicontrol(fg,'Style','Text','Position' ,[310 250 260 016].*WS,...
	'HorizontalAlignment','center', 'String', st.vols{1}.descrip,'FontWeight','bold');
end;


% Positional information
%-----------------------------------------------------------------------
mat = st.vols{1}.premul*st.vols{1}.mat;
Z = spm_imatrix(mat);
Z = Z(7:9);
uicontrol(fg,'Style','Text','Position' ,[310 210 100 016].*WS,...
	'HorizontalAlignment','right', 'String', 'Vox size:');
st.posinf = struct('z',uicontrol(fg,'Style','Text','Position' ,[410 210 160 016].*WS,...
	'HorizontalAlignment','left', 'String', sprintf('%.3g x %.3g x %.3g', Z),'FontWeight','bold'));

O = mat\[0 0 0 1]'; O=O(1:3)';
uicontrol(fg,'Style','Text','Position' ,[310 190 100 016].*WS,...
	'HorizontalAlignment','right', 'String', 'Origin:');
st.posinf.o = uicontrol(fg,'Style','Text','Position' ,[410 190 160 016].*WS,...
	'HorizontalAlignment','left', 'String', sprintf('%.3g %.3g %.3g', O),'FontWeight','bold');

R = spm_imatrix(mat);
R = spm_matrix([0 0 0 R(4:6)]);
R = R(1:3,1:3);

uicontrol(fg,'Style','Text','Position' ,[310 170 100 016].*WS,...
	'HorizontalAlignment','right', 'String', 'Dir Cos:');
tmp2 = sprintf('%+5.3f %+5.3f %+5.3f', R(1,1:3)); tmp2(find(tmp2=='+')) = ' ';
st.posinf.m1 = uicontrol(fg,'Style','Text','Position' ,[410 170 160 016].*WS,...
	'HorizontalAlignment','left', 'String', tmp2,'FontWeight','bold');
tmp2 = sprintf('%+5.3f %+5.3f %+5.3f', R(2,1:3)); tmp2(find(tmp2=='+')) = ' ';
st.posinf.m2 = uicontrol(fg,'Style','Text','Position' ,[410 150 160 016].*WS,...
	'HorizontalAlignment','left', 'String', tmp2,'FontWeight','bold');
tmp2 = sprintf('%+5.3f %+5.3f %+5.3f', R(3,1:3)); tmp2(find(tmp2=='+')) = ' ';
st.posinf.m3 = uicontrol(fg,'Style','Text','Position' ,[410 130 160 016].*WS,...
	'HorizontalAlignment','left', 'String', tmp2,'FontWeight','bold');

tmp = [[R zeros(3,1)] ; 0 0 0 1]*diag([Z 1])*spm_matrix(-O) - mat;
st.posinf.w = uicontrol(fg,'Style','Text','Position' ,[310 110 260 016].*WS,...
	'HorizontalAlignment','center', 'String', '','FontWeight','bold');
if sum(tmp(:).^2)>1e-8,
	set(st.posinf.w, 'String', 'Warning: shears involved');
end;

% Assorted other buttons.
%-----------------------------------------------------------------------
uicontrol(fg,'Style','Frame','Position',[310 30 270 70].*WS);
st.zoomer = uicontrol(fg,'Style','popupmenu' ,'Position',[315 75 125 20].*WS,...
	'String',str2mat('Full Volume','160x160x160mm','80x80x80mm','40x40x40mm','20x20x20mm','10x10x10mm'),...
	'Callback','spm_image(''zoom_in'')','ToolTipString','zoom in by different amounts');
c = 'if get(gco,''Value'')==1, spm_orthviews(''Space''), else, spm_orthviews(''Space'', 1);end;spm_image(''zoom_in'')';
uicontrol(fg,'Style','popupmenu' ,'Position',[315 55 125 20].*WS,...
	'String',str2mat('World Space','Voxel Space'),...
	'Callback',c,'ToolTipString','display in aquired/world orientation');
c = 'if get(gco,''Value'')==1, spm_orthviews(''Xhairs'',''off''), else, spm_orthviews(''Xhairs'',''on''); end;';
uicontrol(fg,'Style','togglebutton','Position',[450 75 125 20].*WS,...
	'String','Hide Crosshairs','Callback',c,'ToolTipString','show/hide crosshairs');
uicontrol(fg,'Style','popupmenu' ,'Position',[450 55 125 20].*WS,...
	'String',str2mat('NN interp','bilin interp','sinc interp'),...
	'Callback','tmp_ = [0 1 -4];spm_orthviews(''Interp'',tmp_(get(gco,''Value'')))',...
	'Value',2,'ToolTipString','interpolation method for displaying images');
st.win = uicontrol(fg,'Style','popupmenu','Position',[315 35 125 20].*WS,...
	'String',str2mat('Auto Window','Manual Window'),'Callback','spm_image(''window'');','ToolTipString','range of voxel intensities displayed');
% uicontrol(fg,'Style','pushbutton','Position',[315 35 125 20].*WS,...
% 	'String','Window','Callback','spm_image(''window'');','ToolTipString','range of voxel intensities % displayed');
st.blobber = uicontrol(fg,'Style','pushbutton','Position',[450 35 125 20].*WS,...
	'String','Add Blobs','Callback','spm_image(''addblobs'');','ToolTipString','superimpose activations');
end;
return;


function my_reset
spm_orthviews('reset');
spm_figure('Clear','Graphics');
return;

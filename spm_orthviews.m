function varargout = spm_orthviews(action,varargin)
% Display Orthogonal Views of a Normalized Image
% FORMAT H = spm_orthviews('Image',filename[,position])
% filename - name of image to display
% area     - position of image
%            -  area(1) - position x
%            -  area(2) - position y
%            -  area(3) - size x
%            -  area(4) - size y
% H        - handle for ortho sections
% FORMAT spm_orthviews('BB',bb)
% bb       - bounding box
%            [loX loY loZ
%             hiX hiY hiZ]
%
% FORMAT spm_orthviews('Redraw')
% Redraws the images
%
% FORMAT spm_orthviews('Reposition',centre)
% centre   - X, Y & Z coordinates of centre voxel
%
% FORMAT spm_orthviews('Space'[,handle])
% handle   - the view to define the space by
% with no arguments - puts things into mm space
%
% FORMAT spm_orthviews('MaxBB')
% sets the bounding box big enough display the whole of all images
%
% FORMAT spm_orthviews('Resolution',res)
% res      - resolution (mm)
%
% FORMAT spm_orthviews('Delete', handle)
% handle   - image number to delete
%
% FORMAT spm_orthviews('Reset')
% clears the orthogonal views
%
% FORMAT spm_orthviews('Pos')
% returns the co-ordinate of the crosshairs in millimetres in the
% standard space.
%
% FORMAT spm_orthviews('Pos', i)
% returns the voxel co-ordinate of the crosshairs in the image in the
% ith orthogonal section.
%
% FORMAT spm_orthviews('Xhairs','off') OR spm_orthviews('Xhairs')
% disables the cross-hairs on the display.
%
% FORMAT spm_orthviews('Xhairs','on')
% enables the cross-hairs.
%
% FORMAT spm_orthviews('Interp',hld)
% sets the hold value to hld (see spm_slice_vol).
%
% FORMAT spm_orthviews('AddBlobs',handle,XYZ,Z,mat)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations (currently in millimeters)
% Z        - blob voxel intensities
% mat      - matrix from millimeters to voxels of blob.
% This method only adds one set of blobs, and displays them using a
% split colour table.
%
% spm_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations (currently in millimeters)
% Z        - blob voxel intensities
% mat      - matrix from millimeters to voxels of blob.
% colour   - the 3 vector containing the colour that the blobs should be
% Several sets of blobs can be added in this way, and it uses full colour.
% Although it may not be particularly attractive on the screen, the colour
% blobs print well.
%
% spm_orthviews('Register',hReg)
% See spm_XYZreg for more information.
%
% FORMAT spm_orthviews('RemoveBlobs',handle)
% Removes all blobs from the image specified by the handle(s).
%
% spm_orthviews('Register',hReg)
% hReg      - Handle of HandleGraphics object to build registry in.
% See spm_XYZreg for more information.
%
%_______________________________________________________________________
% %W% John Ashburner %E%

global st;

if isempty(st), reset_st; end;

spm('Pointer','watch');

if nargin == 0, action = ''; end;
action = lower(action);

switch lower(action),
case 'image',
	H = specify_image(varargin{1});
	if ~isempty(H)
		if length(varargin)>=2, st.vols{H}.area = varargin{2}; end;
		if isempty(st.bb), st.bb = maxbb; end;
		bbox;
		redraw(H);
	end;
	varargout{1} = H;

case 'bb',
	if length(varargin)> 0 & all(size(varargin{1})==[2 3]), st.bb = varargin{1}; end;
	bbox;
	redraw_all;

case 'redraw',
	redraw_all;
	eval(st.callback);
	if isfield(st,'registry'),
		spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
	end;

case 'reposition',
	if length(varargin)<1, tmp = findcent;
	else, tmp = varargin{1}; end;
	if length(tmp)==3, st.centre = tmp(:); end;
	redraw_all;
	eval(st.callback);
	if isfield(st,'registry'),
		spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
	end;

case 'setcoords',
	st.centre = varargin{1};
	st.centre = st.centre(:);
	redraw_all;
	eval(st.callback);

case 'space',
	if length(varargin)<1,
		st.Space = eye(4);
		st.bb = maxbb;
		redraw_all;
	else,
		space(varargin{1});
		redraw_all;
	end;
	bbox;

case 'maxbb',
	st.bb = maxbb;
	bbox;
	redraw_all;

case 'resolution',
	resolution(varargin{1});
	bbox;
	redraw_all;

case 'window',
	if length(varargin)<2,
		win = 'auto';
	elseif length(varargin{2})==2,
		win = varargin{2};
	end;
	for i=valid_handles(varargin{1}),
		st.vols{i}.window = win;
	end;
	redraw(varargin{1});

case 'delete',
	my_delete(varargin{1});

case 'move',
	move(varargin{1},varargin{2});
	% redraw_all;

case 'reset',
	my_reset;

case 'pos',
	if isempty(varargin),
		H = st.centre(:);
	else,
		H = pos(varargin{1});
	end;
	varargout{1} = H;

case 'interp',
	st.hld = varargin{1};
	redraw_all;

case 'xhairs',
	xhairs(varargin{1});

case 'register',
	register(varargin{1});

case 'addblobs',
	addblobs(varargin{1}, varargin{2},varargin{3},varargin{4});
	redraw(varargin{1});

case 'addcolouredblobs',
	addcolouredblobs(varargin{1}, varargin{2},varargin{3},varargin{4},varargin{5});

case 'addimage',
	addimage(varargin{1}, varargin{2});
	redraw(varargin{1});

case 'addcolouredimage',
	addcolouredimage(varargin{1}, varargin{2},varargin{3});

case 'rmblobs',
	rmblobs(varargin{1});
	redraw(varargin{1});

otherwise,
	warning('Unknown action string')
end;

spm('Pointer');
return;


%_______________________________________________________________________
%_______________________________________________________________________
function addblobs(handle, xyz, t, mat)
global st
for i=valid_handles(handle),
	if ~isempty(xyz),
		rcp      = round(xyz);
		dim      = max(rcp,[],2)';
		off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
		vol      = zeros(dim);
		vol(off) = t.*(t > 0);
		vol      = reshape(vol,dim);
		st.vols{i}.blobs=cell(1,1);
		if st.mode == 0,
			axpos = get(st.vols{i}.ax{2}.ax,'Position');
		else,
			axpos = get(st.vols{i}.ax{1}.ax,'Position');
		end;
		ax = axes('Parent',st.fig,'Position',[(axpos(1)+axpos(3)+0.05) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
			'Box','on');
		mx = max([eps max(t)]);
		image([0 mx/32],[0 mx],[1:64]' + 64,'Parent',ax);
		set(ax,'YDir','normal','XTickLabel',[]);
		st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'cbar',ax,'max',mx);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addimage(handle, fname)
global st
for i=valid_handles(handle),
	vol = spm_vol(fname);
	mat = vol.mat;
	st.vols{i}.blobs=cell(1,1);
	if st.mode == 0,
		axpos = get(st.vols{i}.ax{2}.ax,'Position');
	else,
		axpos = get(st.vols{i}.ax{1}.ax,'Position');
	end;
	ax = axes('Parent',st.fig,'Position',[(axpos(1)+axpos(3)+0.05) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
		'Box','on');
	mx = max([eps maxval(vol)]);
	image([0 mx/32],[0 mx],[1:64]' + 64,'Parent',ax);
	set(ax,'YDir','normal','XTickLabel',[]);
	st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'cbar',ax,'max',mx);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredblobs(handle, xyz, t, mat,colour)
global st
for i=valid_handles(handle),
	if ~isempty(xyz),
		rcp      = round(xyz);
		dim      = max(rcp,[],2)';
		off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
		vol      = zeros(dim);
		vol(off) = t.*(t > 0);
		vol      = reshape(vol,dim);
		if ~isfield(st.vols{i},'blobs'),
			st.vols{i}.blobs=cell(1,1);
			bset = 1;
		else,
			bset = length(st.vols{i}.blobs)+1;
		end;
		axpos = get(st.vols{i}.ax{2}.ax,'Position');
		mx = max([eps maxval(vol)]);
		st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'colour',colour);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredimage(handle, fname,colour)
global st
for i=valid_handles(handle),

	vol = spm_vol(fname);
	mat = vol.mat;
	if ~isfield(st.vols{i},'blobs'),
		st.vols{i}.blobs=cell(1,1);
		bset = 1;
	else,
		bset = length(st.vols{i}.blobs)+1;
	end;
	axpos = get(st.vols{i}.ax{2}.ax,'Position');
	mx = max([eps maxval(vol)]);
	st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'colour',colour);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function rmblobs(handle)
global st
for i=valid_handles(handle),
	if isfield(st.vols{i},'blobs'),
		for j=1:length(st.vols{i}.blobs),
			if isfield(st.vols{i}.blobs{j},'cbar') & ishandle(st.vols{i}.blobs{j}.cbar),
				delete(st.vols{i}.blobs{j}.cbar);
			end;
		end;
		st.vols{i} = rmfield(st.vols{i},'blobs');
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function register(hreg)
global st
tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st.fig);
h   = valid_handles(1:24);
if ~isempty(h),
	tmp = st.vols{h(1)}.ax{1}.ax;
	st.registry = struct('hReg',hreg,'hMe', tmp);
	spm_XYZreg('Add2Reg',st.registry.hReg,st.registry.hMe, 'spm_orthviews');
else,
	warning('Nothing to register with');
end;
st.centre = spm_XYZreg('GetCoords',st.registry.hReg);
st.centre = st.centre(:);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function xhairs(arg1),
global st
st.xhairs = 0;
opt = 'on';
if ~strcmp(arg1,'on'),
	opt = 'off';
else,
	st.xhairs = 1;
end;
for i=valid_handles(1:24),
	for j=1:3,
		set(st.vols{i}.ax{j}.lx,'Visible',opt); 
		set(st.vols{i}.ax{j}.ly,'Visible',opt);  
	end; 
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = pos(arg1)
global st
H = [];
for arg1=valid_handles(arg1),
	is = inv(st.vols{arg1}.premul*st.vols{arg1}.mat);
	H = is(1:3,1:3)*st.centre(:) + is(1:3,4);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_reset
global st
if ~isempty(st) & isfield(st,'registry') & ishandle(st.registry.hMe),
	delete(st.registry.hMe); st = rmfield(st,'registry');
end;
my_delete(1:24);
reset_st;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_delete(arg1)
global st
for i=valid_handles(arg1),
	kids = get(st.fig,'Children');
	for j=1:3,
		if any(kids == st.vols{i}.ax{j}.ax),
			set(get(st.vols{i}.ax{j}.ax,'Children'),'DeleteFcn','');
			delete(st.vols{i}.ax{j}.ax);
		end;
	end;
	st.vols{i} = [];
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function resolution(arg1)
global st
res      = arg1/mean(svd(st.Space(1:3,1:3)));
Mat      = diag([res res res 1]);
st.Space = st.Space*Mat;
st.bb    = st.bb/res;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function move(handle,pos)
global st
for handle = valid_handles(handle),
	st.vols{handle}.area = pos;
end;
bbox;
% redraw(valid_handles(handle));
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bb = maxbb
global st
mn = [Inf Inf Inf];
mx = -mn;
for i=valid_handles(1:24),
	bb = [[1 1 1];st.vols{i}.dim(1:3)];
	c = [	bb(1,1) bb(1,2) bb(1,3) 1
		bb(1,1) bb(1,2) bb(2,3) 1
		bb(1,1) bb(2,2) bb(1,3) 1
		bb(1,1) bb(2,2) bb(2,3) 1
		bb(2,1) bb(1,2) bb(1,3) 1
		bb(2,1) bb(1,2) bb(2,3) 1
		bb(2,1) bb(2,2) bb(1,3) 1
		bb(2,1) bb(2,2) bb(2,3) 1]';
	tc = st.Space\(st.vols{i}.premul*st.vols{i}.mat)*c;
	tc = tc(1:3,:)';
	mx = max([tc ; mx]);
	mn = min([tc ; mn]);
end;
bb = [mn ; mx];
return;
%_______________________________________________________________________
%_______________________________________________________________________
function space(arg1)
global st
if ~isempty(st.vols{arg1})
	num = arg1;
	Mat = st.vols{num}.premul(1:3,1:3)*st.vols{num}.mat(1:3,1:3);
	Mat = diag([sqrt(sum(Mat.^2)) 1]);
	Space = (st.vols{num}.premul*st.vols{num}.mat)/Mat;
	bb = [1 1 1;st.vols{num}.dim(1:3)];
	bb = [bb [1;1]];
	bb=bb*Mat';
	bb=bb(:,1:3);
	st.Space  = Space;
	st.bb = bb;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = specify_image(arg1, arg2)
global st
H=[];
ok = 1;
eval('V = spm_vol(arg1);','ok=0;');
if ok == 0,
	fprintf('Can not use image "%s"\n', arg1);
	return;
end;

ii = 1;
while ~isempty(st.vols{ii}), ii = ii + 1; end;

DeleteFcn = ['spm_orthviews(''Delete'',' num2str(ii) ');'];
V.ax = cell(3,1);
for i=1:3,
	ax = axes('Visible','off','DrawMode','fast','Parent',st.fig,'DeleteFcn',DeleteFcn,...
		'YDir','normal');
	d  = image(0,'Tag','Transverse','Parent',ax,...
		'DeleteFcn',DeleteFcn);
	set(ax,'Ydir','normal');
	lx = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	ly = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	if ~st.xhairs,
		set(lx,'Visible','off');
		set(ly,'Visible','off');
	end;
	V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';
st.vols{ii} = V;

H = ii;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bbox
global st
Dims = diff(st.bb)';

TD = Dims([1 2])';
CD = Dims([1 3])';
if st.mode == 0, SD = Dims([3 2])'; else, SD = Dims([2 3])'; end;

un    = get(st.fig,'Units');set(st.fig,'Units','Pixels');sz=get(st.fig,'Position');set(st.fig,'Units',un);
sz    = sz(3:4);
sz(2) = sz(2)-40;

for i=valid_handles(1:24),
	area = st.vols{i}.area(:);
	area = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
	if st.mode == 0,
		sx   = area(3)/(Dims(1)+Dims(3))/1.02;
	else,
		sx   = area(3)/(Dims(1)+Dims(2))/1.02;
	end;
	sy   = area(4)/(Dims(2)+Dims(3))/1.02;
	s    = min([sx sy]);

	offy = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
	sky = s*(Dims(2)+Dims(3))*0.02;
	if st.mode == 0,
		offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
		skx = s*(Dims(1)+Dims(3))*0.02;
	else,
		offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
		skx = s*(Dims(1)+Dims(2))*0.02;
	end;

	DeleteFcn = ['spm_orthviews(''Delete'',' num2str(i) ');'];

	% Transverse
	set(st.vols{i}.ax{1}.ax,'Units','pixels', ...
		'Position',[offx offy s*Dims(1) s*Dims(2)],...
		'Units','normalized','Xlim',[1 TD(1)]+0.5,'Ylim',[1 TD(2)]+0.5,...
		'Visible','on','XTick',[],'YTick',[]);

	% Coronal
	set(st.vols{i}.ax{2}.ax,'Units','Pixels',...
		'Position',[offx offy+s*Dims(2)+sky s*Dims(1) s*Dims(3)],...
		'Units','normalized','Xlim',[1 CD(1)]+0.5,'Ylim',[1 CD(2)]+0.5,...
		'Visible','on','XTick',[],'YTick',[]);

	% Saggital
	if st.mode == 0,
		set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
			'Position',[offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)],...
			'Units','normalized','Xlim',[1 SD(1)]+0.5,'Ylim',[1 SD(2)]+0.5,...
			'Visible','on','XTick',[],'YTick',[]);
	else,
		set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
			'Position',[offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)],...
			'Units','normalized','Xlim',[1 SD(1)]+0.5,'Ylim',[1 SD(2)]+0.5,...
			'Visible','on','XTick',[],'YTick',[]);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function redraw_all
global st
redraw(1:24);
return;
%_______________________________________________________________________
function mx = maxval(vol)
if isstruct(vol),
	mx = -Inf;
	for i=1:vol.dim(3),
		tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
		mx  = max(mx,max(tmp(:)));
	end;
else,
	mx = max(vol(:));
end;

%_______________________________________________________________________
function redraw(arg1)
global st
bb   = st.bb;
Dims = diff(bb)';
is   = inv(st.Space);
cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);

for i = valid_handles(arg1),
	M = st.vols{i}.premul*st.vols{i}.mat;
	TM0 = [	1 0 0 -bb(1,1)
		0 1 0 -bb(1,2)
		0 0 1 -cent(3)
		0 0 0 1];
	TM = inv(TM0*(st.Space\M));
	TD = Dims([1 2]);

	CM0 = [	1 0 0 -bb(1,1)
		0 0 1 -bb(1,3)
		0 1 0 -cent(2)
		0 0 0 1];
	CM = inv(CM0*(st.Space\M));
	CD = Dims([1 3]);

	if st.mode ==0,
		SM0 = [	0 0 1 -bb(1,3)
			0 1 0 -bb(1,2)
			1 0 0 -cent(1)
			0 0 0 1];
		SM = inv(SM0*(st.Space\M)); SD = Dims([3 2]);
	else,
		SM0 = [	0  1 0 -bb(1,2)
			0  0 1 -bb(1,3)
			1  0 0 -cent(1)
			0  0 0 1];
		SM0 = [	0 -1 0 +bb(2,2)
			0  0 1 -bb(1,3)
			1  0 0 -cent(1)
			0  0 0 1];
		SM = inv(SM0*(st.Space\M));
		SD = Dims([2 3]);
	end;

	ok=1;
	eval('imgt  = (spm_slice_vol(st.vols{i},TM,TD,st.hld))'';','ok=0;');
	eval('imgc  = (spm_slice_vol(st.vols{i},CM,CD,st.hld))'';','ok=0;');
	eval('imgs  = (spm_slice_vol(st.vols{i},SM,SD,st.hld))'';','ok=0;');
	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);
	else,
		if strcmp(st.vols{i}.window,'auto'),
			mx = max([max(max(imgt)) max(max(imgc)) max(max(imgs))]);
			mn = min([min(min(imgt)) min(min(imgc)) min(min(imgs))]);
		else,
			mx = st.vols{i}.window(2);
			mn = st.vols{i}.window(1);
			r=min([mn mx]);imgt = max(imgt,r); r=max([mn mx]);imgt = min(imgt,r);
			r=min([mn mx]);imgc = max(imgc,r); r=max([mn mx]);imgc = min(imgc,r);
			r=min([mn mx]);imgs = max(imgs,r); r=max([mn mx]);imgs = min(imgs,r);
		end;

		if isfield(st.vols{i},'blobs'),
			if ~isfield(st.vols{i}.blobs{1},'colour'),
				% Add blobs for display using the split colourmap
				scal = 64/(mx-mn);
				dcoff = -mn*scal;
				imgt = imgt*scal+dcoff;
				imgc = imgc*scal+dcoff;
				imgs = imgs*scal+dcoff;

				vol = st.vols{i}.blobs{1}.vol;
				mat = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				if isfield(st.vols{i}.blobs{1},'max'),
					mx = st.vols{i}.blobs{1}.max;
				else,
					mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.max = mx;
				end;
				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,0)';
				tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,0)';
				tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,0)';

				msk  = find(tmpt>0); imgt(msk) = 64+tmpt(msk)*(64/mx);
				msk  = find(tmpc>0); imgc(msk) = 64+tmpc(msk)*(64/mx);
				msk  = find(tmps>0); imgs(msk) = 64+tmps(msk)*(64/mx);
				colormap([gray(64); pink(64)]);
			else,
				% Add full colour blobs - several sets at once
				scal = 1/(mx-mn);
				dcoff = -mn*scal;
				imgt = repmat(imgt*scal+dcoff,[1,1,3]);
				imgc = repmat(imgc*scal+dcoff,[1,1,3]);
				imgs = repmat(imgs*scal+dcoff,[1,1,3]);

				for j=1:length(st.vols{i}.blobs),
					vol = st.vols{i}.blobs{j}.vol;
					mat = st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
					if isfield(st.vols{i}.blobs{j},'colour'),
						colour = st.vols{i}.blobs{j}.colour;
					else,
						colour = [1 0 0];
					end;
					if isfield(st.vols{i}.blobs{j},'max'),
						mx = st.vols{i}.blobs{j}.max;
					else,
						mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.max = mx;
					end;
					vol  = st.vols{i}.blobs{j}.vol;
					M    = st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
					tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,1)';
					tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,1)';
					tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,1)';

					tmp  = cat(3,tmpt*(colour(1)/mx),tmpt*(colour(2)/mx),tmpt*(colour(3)/mx));
					imgt = (repmat((1-tmpt/mx),[1 1 3]).*imgt+tmp);
					tmp = find(imgt<0); imgt(tmp)=0; tmp = find(imgt>1); imgt(tmp)=1;

					tmp  = cat(3,tmpc*(colour(1)/mx),tmpc*(colour(2)/mx),tmpc*(colour(3)/mx));
					imgc = (repmat((1-tmpc/mx),[1 1 3]).*imgc+tmp);
					tmp = find(imgc<0); imgc(tmp)=0; tmp = find(imgc>1); imgc(tmp)=1;

					tmp  = cat(3,tmps*(colour(1)/mx),tmps*(colour(2)/mx),tmps*(colour(3)/mx));
					imgs = (repmat((1-tmps/mx),[1 1 3]).*imgs+tmp);
					tmp = find(imgs<0); imgs(tmp)=0; tmp = find(imgs>1); imgs(tmp)=1;
				end;
			end;
		else,
			scal = 64/(mx-mn);
			dcoff = -mn*scal;
			imgt = imgt*scal+dcoff;
			imgc = imgc*scal+dcoff;
			imgs = imgs*scal+dcoff;
		end;

		callback = 'spm_orthviews(''Reposition'');';

		set(st.vols{i}.ax{1}.d,'ButtonDownFcn',callback, 'Cdata',imgt);
		set(st.vols{i}.ax{1}.lx,'ButtonDownFcn',callback,...
			'Xdata',[0 TD(1)],'Ydata',[1 1]*(cent(2)-bb(1,2)));
		set(st.vols{i}.ax{1}.ly,'ButtonDownFcn',callback,...
			'Ydata',[0 TD(2)],'Xdata',[1 1]*(cent(1)-bb(1,1)));

		set(st.vols{i}.ax{2}.d,'ButtonDownFcn',callback, 'Cdata',imgc);
		set(st.vols{i}.ax{2}.lx,'ButtonDownFcn',callback,...
			'Xdata',[0 CD(1)],'Ydata',[1 1]*(cent(3)-bb(1,3)));
		set(st.vols{i}.ax{2}.ly,'ButtonDownFcn',callback,...
			'Ydata',[0 CD(2)],'Xdata',[1 1]*(cent(1)-bb(1,1)));

		set(st.vols{i}.ax{3}.d,'ButtonDownFcn',callback,'Cdata',imgs);
		if st.mode ==0,
			set(st.vols{i}.ax{3}.lx,'ButtonDownFcn',callback,...
				'Xdata',[0 SD(1)],'Ydata',[1 1]*(cent(2)-bb(1,2)));
			set(st.vols{i}.ax{3}.ly,'ButtonDownFcn',callback,...
				'Ydata',[0 SD(2)],'Xdata',[1 1]*(cent(3)-bb(1,3)));
		else,
			set(st.vols{i}.ax{3}.lx,'ButtonDownFcn',callback,...
				'Xdata',[0 SD(1)],'Ydata',[1 1]*(cent(3)-bb(1,3)));
			set(st.vols{i}.ax{3}.ly,'ButtonDownFcn',callback,...
				'Ydata',[0 SD(2)],'Xdata',[1 1]*(bb(2,2)-cent(2)));
		end;
	end;
end;
drawnow;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function centre = findcent
global st
obj    = get(st.fig,'CurrentObject');
centre = [];
cent   = [];
cp     = [];
for i=valid_handles(1:24),
	for j=1:3,
		if ~isempty(obj),
			if any([st.vols{i}.ax{j}.d  ...
				st.vols{i}.ax{j}.lx ...
				st.vols{i}.ax{j}.ly]== obj)
				cp = get(get(obj,'Parent'),'CurrentPoint');
			elseif (st.vols{i}.ax{j}.ax == obj),
				cp = get(obj,'CurrentPoint');
			end;
		end;
		if ~isempty(cp),
			cp   = cp(1,1:2);
			is   = inv(st.Space);
			cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
			switch j,
				case 1,
				cent([1 2])=[cp(1)+st.bb(1,1) cp(2)+st.bb(1,2)];
				case 2,
				cent([1 3])=[cp(1)+st.bb(1,1) cp(2)+st.bb(1,3)];
				case 3,
				if st.mode ==0,
					cent([3 2])=[cp(1)+st.bb(1,3) cp(2)+st.bb(1,2)];
				else,
					cent([2 3])=[st.bb(2,2)-cp(1) cp(2)+st.bb(1,3)];
				end;
			end;
			break;
		end;
	end;
	if ~isempty(cent), break; end;
end;
if ~isempty(cent), centre = st.Space(1:3,1:3)*cent(:) + st.Space(1:3,4); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function handles = valid_handles(handles)
global st;
handles = handles(:)';
handles = handles(find(handles<=24 & handles>=1 & ~rem(handles,1)));
for h=handles,
	if isempty(st.vols{h}), handles(find(handles==h))=[]; end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function reset_st
global st
fig     = spm_figure('FindWin','Graphics');
bb      = [ [-78 78]' [-112 76]' [-50 85]' ];
st      = struct('n', 0, 'vols',[], 'bb',bb,'Space',eye(4),'centre',[0 0 0],'callback',';','xhairs',1,'hld',1,'fig',fig,'mode',1);
st.vols = cell(24,1);
return;

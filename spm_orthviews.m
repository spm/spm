function H = spm_orthviews(action,arg1,arg2,arg3)
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
% FORMAT spm_orthviews('Resolution',fac)
% fac      - increase the voxel sizes by factor fac
%
% FORMAT spm_orthviews('Delete', handle)
% handle   - image number to delete
%
% FORMAT spm_orthviews('Reset')
% clears the orthogonal views
%_______________________________________________________________________
% %W% John Ashburner %E%


global centre;

if (nargin==0)
	spm_orthviews('BB');
	return;
end

global st;
fig = spm_figure('FindWin','Graphics');

if isempty(st)
	bb     = [ [-78 78]' [-112 76]' [-50 85]' ];
	st = struct('n', 0, 'vols',[], 'bb',bb,'Space',eye(4),'centre',[0 0 0],'callback','disp(centre)');
	st.vols = cell(24,1);
end

action = lower(action);

if strcmp(action,'image')
	if (nargin<2), return; end;

	ok = 1;
	eval('[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(arg1);','ok=0;');
	eval('M = spm_get_space(arg1);','ok=0;');
	if ok == 0,
		fprintf('Can not use image "%s"\n', arg1);
		return;
	end

	D = [DIM VOX SCALE TYPE OFFSET];

	if nargin>2,
		area = arg2;
	else
		area = [0. 0. 1. 1.];
	end

	ii = 1;
	while ~isempty(st.vols{ii})
		ii = ii + 1;
	end

	st.vols{ii} = struct('V',0,'D',D,'M',M,'name',arg1,...
		'area',area,'ax',[]);
	ax = cell(3,1);
	DeleteFcn = ['spm_orthviews(''Delete'',' num2str(ii) ');'];
	for i=1:3
		ax = axes('Visible','off','DrawMode','fast','Parent',fig,'DeleteFcn',DeleteFcn,...
			'YDir','normal');
		d  = image(0,'Tag','Transverse','Parent',ax,...
			'DeleteFcn',DeleteFcn);
		set(ax,'Ydir','normal');
		lx = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
		ly = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
		axx{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
	end
	st.vols{ii}.ax = axx;

	H = ii;
	if isempty(st.bb),
		spm_orthviews('Space', H);
	else
		spm_orthviews('BB',st.bb);
	end
end

if strcmp(action,'bb')
	if nargin == 1,
		bb = st.bb;
	else
		bb = arg1;
	end;


	Dims = diff(bb)';

	TD = Dims([1 2])';
	CD = Dims([1 3])';
	SD = Dims([3 2])';

	un=get(fig,'Units');set(fig,'Units','Pixels');sz=get(fig,'Position');set(fig,'Units',un);
	sz = sz(3:4);
	sz(2) = sz(2)-40;

	for i=1:24
		if ~isempty(st.vols{i})

			area = st.vols{i}.area(:);
			area = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
			sx = area(3)/(Dims(1)+Dims(3));
			sy = area(4)/(Dims(2)+Dims(3));
			s = min([sx sy]);
			offx = (area(3)-(Dims(1)+Dims(3))*s)/2 + area(1);
			offy = (area(4)-(Dims(2)+Dims(3))*s)/2 + area(2);

			DeleteFcn = ['spm_orthviews(''Delete'',' num2str(i) ');'];

			% Transverse
			set(st.vols{i}.ax{1}.ax,'Units','pixels', ...
				'Position',[offx offy s*Dims(1) s*Dims(2)],...
				'Units','normalized','Xlim',[1 TD(1)]+0.5,'Ylim',[1 TD(2)]+0.5,...
				'Visible','off');

			% Coronal
			set(st.vols{i}.ax{2}.ax,'Units','Pixels',...
				'Position',[offx offy+s*Dims(2) s*Dims(1) s*Dims(3)],...
				'Units','normalized','Xlim',[1 CD(1)]+0.5,'Ylim',[1 CD(2)]+0.5,...
				'Visible','off');

			% Saggital
			set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
				'Position',[offx+s*Dims(1) offy s*Dims(3) s*Dims(2)],...
				'Units','normalized','Xlim',[1 SD(1)]+0.5,'Ylim',[1 SD(2)]+0.5,...
				'Visible','off');

		end
	end
	st.bb = bb;
	spm_orthviews('Reposition',st.centre);
end


if strcmp(action,'reposition')
	centre = [];
	if (nargin == 1)
		cp = [];
		obj = get(fig,'CurrentObject');
		a = 0;
		for i=1:24
			if ~isempty(st.vols{i})
				for j=1:3
					if any([st.vols{i}.ax{j}.d  ...
						st.vols{i}.ax{j}.lx ...
						st.vols{i}.ax{j}.ly]== obj)
						cp = get(get(obj,'Parent'),'CurrentPoint');
					elseif (st.vols{i}.ax{j}.ax == obj)
						cp = get(obj,'CurrentPoint');
					end
					if ~isempty(cp)
						cp = cp(1,1:2);
						centre = st.centre;
						switch j
							case 1,
							centre([1 2])=[cp(1)+st.bb(1,1) cp(2)+st.bb(1,2)];
							case 2,
							centre([1 3])=[cp(1)+st.bb(1,1) cp(2)+st.bb(1,3)];
							case 3,
							centre([3 2])=[cp(1)+st.bb(1,3) cp(2)+st.bb(1,2)];
						end
						break;
					end
				end
				if ~isempty(centre), break; end;
			end
		end
		if isempty(centre), return; end;
	else
		centre = arg1;
	end

	bb = st.bb;
	Dims = diff(bb)';

	for i=1:24
		if ~isempty(st.vols{i})
			M=st.vols{i}.M;
			TM0 = [
				1 0 0 -bb(1,1)
				0 1 0 -bb(1,2)
				0 0 1 -centre(3)
				0 0 0 1];

			CM0 = [
				1 0 0 -bb(1,1)
				0 0 1 -bb(1,3)
				0 1 0 -centre(2)
				0 0 0 1];
	
			SM0 = [
				0 0 1 -bb(1,3)
				0 1 0 -bb(1,2)
				1 0 0 -centre(1)
				0 0 0 1];

			TM = inv(TM0*(st.Space\M)); TD = Dims([1 2]);
			CM = inv(CM0*(st.Space\M)); CD = Dims([1 3]);
			SM = inv(SM0*(st.Space\M)); SD = Dims([3 2]);

			mapped = 1;
			if isempty(st.vols{i}.V) | st.vols{i}.V == 0
				mapped = 1;
				eval('V = spm_map_vol(st.vols{i}.name,st.vols{i}.D);','mapped = 0;');
				if isempty(st.vols{i}.V),
					st.vols{i}.V = V;
				end
			else
				V = st.vols{i}.V;
			end
			if (mapped == 0)
				fprintf('Image "%s" can not be resampled\n', st.vols{i}.name);
			else
				imgt  = (spm_slice_vol(V,TM,TD,1))';
				imgc  = (spm_slice_vol(V,CM,CD,1))';
				imgs  = (spm_slice_vol(V,SM,SD,1))';
	
	
				if ~isempty(st.vols{i}.V) & st.vols{i}.V==0, spm_unmap_vol(V); end;
	
				scal = 64/max([max(max(imgt)) max(max(imgc)) max(max(imgs))]);
	
				posn = [centre(1)+bb(1,1) centre(2)+bb(1,2) centre(3)+bb(1,3) 1]';
				posn = [centre(1)-bb(1,1) centre(2)-bb(1,2) centre(3)-bb(1,3) 1]';
				posn = posn(1:3);
	
				callback = 'spm_orthviews(''Reposition'');';

				set(st.vols{i}.ax{1}.d,'ButtonDownFcn',callback, 'Cdata',imgt*scal);
				set(st.vols{i}.ax{1}.lx,'ButtonDownFcn',callback,...
					'Xdata',[0 TD(1)],'Ydata',[1 1]*posn(2));
				set(st.vols{i}.ax{1}.ly,'ButtonDownFcn',callback,...
					'Ydata',[0 TD(2)],'Xdata',[1 1]*posn(1));

				set(st.vols{i}.ax{2}.d,'ButtonDownFcn',callback, 'Cdata',imgc*scal);
				set(st.vols{i}.ax{2}.lx,'ButtonDownFcn',callback,...
					'Xdata',[0 CD(1)],'Ydata',[1 1]*posn(3));
				set(st.vols{i}.ax{2}.ly,'ButtonDownFcn',callback,...
					'Ydata',[0 CD(2)],'Xdata',[1 1]*posn(1));

				set(st.vols{i}.ax{3}.d,'ButtonDownFcn',callback,'Cdata',imgs*scal);
				set(st.vols{i}.ax{3}.lx,'ButtonDownFcn',callback,...
					'Xdata',[0 SD(1)],'Ydata',[1 1]*posn(2));
				set(st.vols{i}.ax{3}.ly,'ButtonDownFcn',callback,...
					'Ydata',[0 SD(2)],'Xdata',[1 1]*posn(3));
				drawnow;
			end
		end
	end
	st.centre = centre;
	eval(st.callback);
end



if (strcmp(action,'space'))
	Space = eye(4);
	bb = [ [-64 64]' [-104 68]' [-28 72]' ];
	
	if (nargin>1)
		if ~isempty(st.vols{arg1})
			num = arg1;
			Mat = st.vols{num}.M(1:3,1:3);
			Mat = diag([sqrt(sum(Mat.^2)) 1]);
			Space = st.vols{num}.M/Mat;
			bb = [1 1 1;st.vols{num}.D(1:3)];
			bb = [bb [1;1]];
			bb=bb*Mat';
			bb=bb(:,1:3);
		end
	end
	st.centre = (Space\st.Space)*[st.centre';1];
	st.centre = st.centre(1:3)';
	st.Space  = Space;
	st.bb = bb;
	spm_orthviews('BB',bb);
end

if (strcmp(action,'maxbb'))
	mn = [Inf Inf Inf];
	mx = -mn;
	for i=1:24
		if ~isempty(st.vols{i})
			bb = [[1 1 1];st.vols{i}.D(1:3)];
			c = [	bb(1,1) bb(1,2) bb(1,3) 1
				bb(1,1) bb(1,2) bb(2,3) 1
				bb(1,1) bb(2,2) bb(1,3) 1
				bb(1,1) bb(2,2) bb(2,3) 1
				bb(2,1) bb(1,2) bb(1,3) 1
				bb(2,1) bb(1,2) bb(2,3) 1
				bb(2,1) bb(2,2) bb(1,3) 1
				bb(2,1) bb(2,2) bb(2,3) 1]';
			tc = st.Space*st.vols{i}.M*c;
			tc = tc(1:3,:)';
			mx = max([tc ; mx]);
			mn = min([tc ; mn]);
		end
	end
	spm_orthviews('BB',[mn ; mx]);
end


if (strcmp(action,'resolution'))
	if nargin ~= 2, return; end
	res = arg1;
	Mat = diag([res res res 1]);
	st.Space = st.Space*Mat;
	st.bb = st.bb/res;
	st.centre = st.centre/res;
	spm_orthviews('BB',st.bb);
end

if (strcmp(action,'delete'))
	if nargin ~= 2, return; end
	arg1 = arg1(find(arg1>=1 & arg1<=24));
	for i=arg1
		if ~isempty(st.vols{i})
			kids = get(fig,'Children');
			for j=1:3,
				if any(kids == st.vols{i}.ax{j}.ax),
					set(get(st.vols{i}.ax{j}.ax,'Children'),'DeleteFcn','');
					delete(st.vols{i}.ax{j}.ax);
				end
			end
			if ~isempty(st.vols{i}.V) & all(st.vols{i}.V~=0), spm_unmap(st.vols{i}.V); end;
			st.vols{i} = [];
		end
	end
end

if (strcmp(action,'replace'))
	if nargin ~= 3, return; end
	image = arg1;
	del   = arg2;
	area = st.vols{del}.area;
	spm_orthviews('Delete',del);
	spm_orthviews('Image',image,area);
end

if (strcmp(action,'reset'))
	if nargin ~= 1, return; end
	spm_orthviews('Delete',1:24);
	st = [];
end


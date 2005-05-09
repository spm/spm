function ret = spm_ov_roi(varargin)

% ROI tool - plugin for spm_orthviews
%
% With ROI tool it is possible to create new or modify existing mask images
% interactively. ROI tool can be launched via the spm_orthviews image
% context menu.
% While ROI tool is active, mouse buttons have the following functions:
% left       Reposition crosshairs
% middle     Perform ROI tool operation according to selected edit mode at
%            crosshair position
% right      context menu
% Menu options and prompts explained:
% Launch     Initialise ROI tool in current image
%            'Load existing ROI image? (yes/no)' 
%              If you want to modify an existing mask image (e.g. mask.img from
%              a fMRI analysis), press 'yes'. You will then be prompted to
%            'Select ROI image'
%              This is the image that will be loaded as initial ROI.
%              If you want to create a new ROI image, you will first be
%              prompted to
%            'Select image defining ROI space'
%              The image dimensions, voxel sizes and slice orientation will
%              be read from this image. Thus you can edit a ROI based on a
%              image with a resolution and slice orientation different from
%              the underlying displayed image.
%            'ROI filename'
%              Defaults to 'PWD/roitool.img'. This is the filename for your
%              ROI image when you select the save menu option.
% Edit mode  Operation performed when pressing the middle mouse button.
%            'rect' 
%              You will be prompted for the size of a box in voxels. The box
%              will be placed centered at crosshair position. 
%            'toggle'
%              Toggle the selection state of the voxel at crosshair position. 
%            'unset'
%              Unset the voxel at crosshair position.
% Threshold  You will be prompted to enter a [min max] threshold. Only
%            those voxels in the ROI image where the intensities of the
%            underlying image are within the [min max] range will survive
%            this operation.            
% Clear      Clear ROI, but keep ROI space information
% Save       Save ROI image
% Quit       Quit ROI tool

% Note: This plugin depends on the blobs set by spm_orthviews('addblobs',...) 
% They should not be removed while ROI tool is active and no other blobs be
% added. This restriction may be removed when switching to MATLAB 6.x and
% using the 'alpha' property to overlay blobs onto images.

global st;
if isempty(st)
  error('roi: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
  error('roi: Wrong number of arguments. Usage: spm_orthviews(''roi'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};

toset = [];
toclear = [];
update_roi = 0;

switch cmd
  case 'init'
    spm('pointer','watch');
    Vroi = spm_vol(varargin{3});
    if nargin > 3
      froi = varargin{4};
    else
      froi = [];
    end;
    if isempty(froi)
      roi = spm_read_vols(Vroi);
      [x y z] = ndgrid(1:Vroi.dim(1),1:Vroi.dim(2),1:Vroi.dim(3));
      xyz = [x(:)'; y(:)'; z(:)'];
    else
      roi = zeros(Vroi.dim(1:3));
      Vroi.fname = froi;
      Vroi.dim(4) = 4;
      xyz = [];
    end;

    clear x y z

    % draw a frame only if ROI volume different from underlying GM volume
    if any(Vroi.dim(1:3)-st.vols{volhandle}.dim(1:3))| ...
	  any(Vroi.mat(:)-st.vols{volhandle}.mat(:))
      [xx1 yx1 zx1] = ndgrid(1            , 1:Vroi.dim(2), 1:Vroi.dim(3));
      [xx2 yx2 zx2] = ndgrid(Vroi.dim(1)  , 1:Vroi.dim(2), 1:Vroi.dim(3));
      [xy1 yy1 zy1] = ndgrid(1:Vroi.dim(1), 1            , 1:Vroi.dim(3));
      [xy2 yy2 zy2] = ndgrid(1:Vroi.dim(1), Vroi.dim(2)  , 1:Vroi.dim(3));
      [xz1 yz1 zz1] = ndgrid(1:Vroi.dim(1), 1:Vroi.dim(2), 1);
      [xz2 yz2 zz2] = ndgrid(1:Vroi.dim(1), 1:Vroi.dim(2), Vroi.dim(3));
      
      fxyz = [xx1(:)' xx2(:)' xy1(:)' xy2(:)' xz1(:)' xz2(:)'; ...
              yx1(:)' yx2(:)' yy1(:)' yy2(:)' yz1(:)' yz2(:)'; ...
	      zx1(:)' zx2(:)' zy1(:)' zy2(:)' zz1(:)' zz2(:)'];
      clear xx1 yx1 zx1 xx2 yx2 zx2 xy1 yy1 zy1 xy2 yy2 zy2 xz1 yz1 zz1 xz2 yz2 zz2
    else
      fxyz  = [];
    end;

    for k=1:3
      cb{k}=get(st.vols{volhandle}.ax{k}.ax,'ButtonDownFcn');
      set(st.vols{volhandle}.ax{k}.ax,'ButtonDownFcn',...
	  ['switch get(gcf,''SelectionType'')',...
	    'case ''normal'', spm_orthviews(''Reposition'');',...
	    'case ''extend'', spm_orthviews(''roi'',''edit'',', ...
	      num2str(volhandle), ');',...
	    'case ''alt'', spm_orthviews(''context_menu'',''ts'',1);',...
	    'end;']);
    end;

    st.vols{volhandle}.roi = struct('Vroi', Vroi, 'xyz', xyz, 'roi', ...
	roi, 'hroi', 1, 'fxyz', fxyz, ...
	'hframe', 1, 'mode', 'rect', 'thresh', [60 140], 'rect', ...
	[4 4 4], 'cb', []);
    st.vols{volhandle}.roi.cb = cb;

    if ~isempty(st.vols{volhandle}.roi.fxyz)
      if isfield(st.vols{volhandle}, 'blobs')
        st.vols{volhandle}.roi.hframe = prod(size(st.vols{volhandle}.blobs))+1;
      end;
      spm_orthviews('addcolouredblobs', volhandle, ...
	  st.vols{volhandle}.roi.fxyz, ones(size(st.vols{volhandle}.roi.fxyz,2),1), ... 
	  st.vols{volhandle}.roi.Vroi.mat,[1 .5 .5]);
      st.vols{volhandle}.blobs{st.vols{volhandle}.roi.hframe}.max=1.3;
    end;
    update_roi=1;

  case 'edit'
    spm('pointer','watch');
    pos = round(inv(st.vols{volhandle}.roi.Vroi.mat)* ...
	[spm_orthviews('pos'); 1]); 
    switch st.vols{volhandle}.roi.mode
      case 'toggle'
	if all(pos(1:3)>0) & all(pos(1:3)<=st.vols{volhandle}.roi.Vroi.dim(1:3)')
	  if st.vols{volhandle}.roi.roi(pos(1),pos(2),pos(3)) 
	    toclear = pos(1:3);
          else
            toset = pos(1:3);
          end;
        end;
  
      case 'rect'    
	tmp = round((st.vols{volhandle}.roi.rect-1)/2);
	[sx sy sz] = meshgrid(-tmp(1):tmp(1), -tmp(2):tmp(2), -tmp(3):tmp(3));
	sel = [sx(:)';sy(:)';sz(:)']+repmat(pos(1:3), 1,prod(2*tmp+1));
	toset = sel(:, (all(sel>0) &...
	    sel(1,:)<=st.vols{volhandle}.roi.Vroi.dim(1) & ...
	    sel(2,:)<=st.vols{volhandle}.roi.Vroi.dim(2) & ...
	    sel(3,:)<=st.vols{volhandle}.roi.Vroi.dim(3)));
	
      case 'unset'
	if all(pos(1:3)>0) & all(pos(1:3)<=st.vols{volhandle}.roi.Vroi.dim(1:3)')
          toclear = pos(1:3);
        end;

      otherwise
	warning('spm_orthviews(''roi'',''edit'', volhandle):Unknown edit mode');
    end;
    update_roi = 1;
  
  case 'thresh'
    spm('pointer','watch');
    rind = find(st.vols{volhandle}.roi.roi);
    [x y z]=ind2sub(st.vols{volhandle}.roi.Vroi.dim(1:3),rind);
    tmp = round(inv(st.vols{volhandle}.mat) * ...
	st.vols{volhandle}.roi.Vroi.mat*[x'; y'; z'; ones(size(x'))]); 
    dat = spm_sample_vol(st.vols{volhandle}, ...
	tmp(1,:), tmp(2,:), tmp(3,:), 0);
    sel = ~((st.vols{volhandle}.roi.thresh(1) < dat) & ...
	(dat < st.vols{volhandle}.roi.thresh(2)));
    toclear =  [x(sel)'; y(sel)'; z(sel)'];
    update_roi = 1;

  case 'connect'
    spm('pointer','watch');    
    pos = round(inv(st.vols{volhandle}.roi.Vroi.mat)* ...
	[spm_orthviews('pos'); 1]); 
    pos = pos(1:3);
    toclear = st.vols{volhandle}.roi.xyz;
    toset = zeros(size(st.vols{volhandle}.roi.xyz));
    curset = 1;
    vis = zeros(size(st.vols{volhandle}.roi.roi));
    pstck = stack('init',size(st.vols{volhandle}.roi.xyz,2));
    nbhood = [ -1  1  0  0  0  0; ...
                0  0 -1  1  0  0; ...
		0  0  0  0 -1  1];
    if st.vols{volhandle}.roi.roi(pos(1),pos(2),pos(3))
      toset(:,curset) = pos;
      curset = curset+1;
      vis(pos(1),pos(2),pos(3)) = 1;
      for k = 1:size(nbhood,2)
        cpos = pos + nbhood(:,k);
        if st.vols{volhandle}.roi.roi(cpos(1), cpos(2), cpos(3)) & ...
           ~vis(cpos(1), cpos(2), cpos(3)) & all(cpos>0) & all(cpos<=st.vols{volhandle}.dim(1:3)') 
	  pstck=stack('push',pstck,cpos);
          vis(cpos(1), cpos(2), cpos(3)) = 1;
        end;
      end;
      while ~stack('isempty',pstck)
        [pstck pos] = stack('pop', pstck);
        toset(:,curset) = pos;
        curset = curset+1;
        vis(pos(1),pos(2),pos(3)) = 1;
        for k = 1:size(nbhood,2)
          cpos = pos + nbhood(:,k);
          if st.vols{volhandle}.roi.roi(cpos(1), cpos(2), cpos(3)) & ...
             ~vis(cpos(1), cpos(2), cpos(3)) & all(cpos>0) & all(cpos<=st.vols{volhandle}.dim(1:3)')
  	    pstck = stack('push',pstck,cpos);
            vis(cpos(1), cpos(2), cpos(3)) = 1;
          end;
        end;
      end;
    end;
    toset = toset(:,1:curset-1);
    update_roi = 1;

  case 'clear'
    spm('pointer','watch');
    st.vols{volhandle}.roi.roi = zeros(size(st.vols{volhandle}.roi.roi));
    st.vols{volhandle}.roi.xyz=[];
    update_roi = 1;
    
  case 'save'
    spm('pointer','watch');
    spm_write_vol(st.vols{volhandle}.roi.Vroi, st.vols{volhandle}.roi.roi);
  
  case 'set' % no type checking!!
    if isfield(st.vols{volhandle}.roi, varargin{3})
      st.vols{volhandle}.roi=setfield(st.vols{volhandle}.roi, varargin{3}, varargin{4});
    else
      warning(sprintf('spm_ov_roi(''set'', ...): Unknown field %s', ...
	  varargin{3}));
    end;
    
  case 'redraw'
    % do nothing
  
  %-------------------------------------------------------------------------
  % Context menu and callbacks
  case 'context_menu'  
    item0 = uimenu(varargin{3}, 'Label', 'ROI tool');
      item1 = uimenu(item0, 'Label', 'Launch', 'Callback', ...
	  ['feval(''spm_ov_roi'',''context_init'', ', ...
	    num2str(volhandle), ');'], 'Tag', ['ROI_0_', num2str(volhandle)]);
      item2 = uimenu(item0, 'Label', 'Edit mode', ...
	  'Visible', 'off', 'Tag', ['ROI_1_', num2str(volhandle)]);
        item2_1 = uimenu(item2, 'Label', 'rect', 'Callback', ...
	    ['feval(''spm_ov_roi'',''context_edit'',', ...
	      num2str(volhandle), ',''rect'');'], ...
	    'Tag', ['ROI_EDIT_', num2str(volhandle)]);
	item2_2 = uimenu(item2, 'Label', 'toggle', 'Callback', ...
	    ['feval(''spm_ov_roi'',''context_edit'',', ...
	      num2str(volhandle), ',''toggle'');'], ...
	    'Tag', ['ROI_EDIT_', num2str(volhandle)]);
	item2_3 = uimenu(item2, 'Label', 'unset', 'Callback', ...
	    ['feval(''spm_ov_roi'',''context_edit'',', ...
	      num2str(volhandle), ',''unset'');'], ...
	    'Tag', ['ROI_EDIT_', num2str(volhandle)]); 
      item3 = uimenu(item0, 'Label', 'Threshold', 'Callback', ...
	  ['feval(''spm_ov_roi'',''context_thresh'', ', ...
	    num2str(volhandle), ');'], 'Visible', 'off', ...
	  'Tag', ['ROI_1_', num2str(volhandle)]);
      item4 = uimenu(item0, 'Label', 'Select connected', 'Callback', ...
	  ['feval(''spm_ov_roi'',''connect'', ', ...
	    num2str(volhandle), ');'], 'Visible', 'off', ...
	  'Tag', ['ROI_1_', num2str(volhandle)]);
      item5 = uimenu(item0, 'Label', 'Clear', 'Callback', ...
	  ['feval(''spm_ov_roi'',''clear'', ', ...
	    num2str(volhandle), ');'], 'Visible', 'off', ...
	  'Tag', ['ROI_1_', num2str(volhandle)]);
      item6 = uimenu(item0, 'Label', 'Save', 'Callback', ...
	  ['feval(''spm_ov_roi'',''save'', ', ...
	    num2str(volhandle), ');'], 'Visible', 'off', ...
	  'Tag', ['ROI_1_', num2str(volhandle)]);
      item7 = uimenu(item0, 'Label', 'Quit', 'Callback', ...
	  ['feval(''spm_ov_roi'',''context_quit'', ', ...
	    num2str(volhandle), ');'], 'Visible', 'off', ...
	  'Tag', ['ROI_1_', num2str(volhandle)]);
    
  case 'context_init'
    Finter = spm_figure('FindWin', 'Interactive');
    spm_input('!DeleteInputObj',Finter);
    usefile = spm_input('Load existing ROI image?','!+1','b','yes|no',[1 0],1);
    if usefile
      imfname = spm_select(1, 'image', 'Select ROI image');
      roifname = [];
    else
      imfname = spm_select(1, 'image', 'Select image defining ROI space');
      [p n e v] = fileparts(imfname);
      roifname = fullfile(p,['roitool' e v]);
      roifname = spm_input('ROI filename','!+1','s',roifname);
    end;
    feval('spm_ov_roi','init',volhandle,imfname,roifname);
    obj = findobj(0, 'Tag',  ['ROI_1_', num2str(volhandle)]);
    set(obj, 'Visible', 'on');
    obj = findobj(0, 'Tag',  ['ROI_0_', num2str(volhandle)]);
    set(obj, 'Visible', 'off');
    obj = findobj(0, 'Tag', ['ROI_EDIT_', num2str(volhandle)]);
    set(obj, 'Checked', 'off');
    obj = findobj(0, 'Tag', ['ROI_EDIT_', num2str(volhandle)], 'Label', st.vols{volhandle}.roi.mode);
    set(obj, 'Checked', 'on');
    spm_input('!DeleteInputObj',Finter);
    
  case 'context_edit'
    feval('spm_ov_roi', 'set', volhandle, 'mode', varargin{3});
    obj = findobj(0, 'Tag', ['ROI_EDIT_', num2str(volhandle)]);
    set(obj, 'Checked', 'off');
    obj = findobj(0, 'Tag', ['ROI_EDIT_', num2str(volhandle)], 'Label', st.vols{volhandle}.roi.mode);
    set(obj, 'Checked', 'on');
    if strcmp(varargin{3}, 'rect')
      Finter = spm_figure('FindWin', 'Interactive');
      spm_input('!DeleteInputObj',Finter);
      rect = spm_input('Selection size {vx vy vz}','!+1','e', ...
	  num2str(st.vols{volhandle}.roi.rect), [1 3]);
      spm_input('!DeleteInputObj',Finter);
      st.vols{volhandle}.roi.rect = rect;
    end;
    
  case 'context_thresh'
    Finter = spm_figure('FindWin', 'Interactive');
    spm_input('!DeleteInputObj',Finter);
    thresh = spm_input('Threshold  {min max}','!+1','e', ...
	num2str(st.vols{volhandle}.roi.thresh), [1 2]);
    spm_input('!DeleteInputObj',Finter);
    st.vols{volhandle}.roi.thresh = thresh;
    feval('spm_ov_roi', 'thresh', volhandle);
    
  case 'context_quit'
    obj = findobj(0, 'Tag',  ['ROI_1_', num2str(volhandle)]);
    set(obj, 'Visible', 'off');
    obj = findobj(0, 'Tag',  ['ROI_0_', num2str(volhandle)]);
    set(obj, 'Visible', 'on');
    spm_orthviews('rmblobs', volhandle);
    for k=1:3
      set(st.vols{volhandle}.ax{k}.ax,'ButtonDownFcn', st.vols{volhandle}.roi.cb{k});
    end;
    st.vols{volhandle} = rmfield(st.vols{volhandle}, 'roi');

  otherwise    
    fprintf('spm_orthviews(''roi'', ...): Unknown action %s', cmd);
  end;

if update_roi
  % clear first, then set (needed for connect operation)
  if ~isempty(toclear)
    itoclear = sub2ind(st.vols{volhandle}.roi.Vroi.dim(1:3), ...
	      toclear(1,:), toclear(2,:), toclear(3,:)); 
    st.vols{volhandle}.roi.roi(itoclear) = 0;
    if ~isempty(st.vols{volhandle}.roi.xyz)
      st.vols{volhandle}.roi.xyz = setdiff(st.vols{volhandle}.roi.xyz',toclear','rows')';
    else
      st.vols{volhandle}.roi.xyz = [];
    end;
  end;

  if ~isempty(toset)
    itoset = sub2ind(st.vols{volhandle}.roi.Vroi.dim(1:3), ...
            toset(1,:), toset(2,:), toset(3,:)); 
    st.vols{volhandle}.roi.roi(itoset) = 1;
    if ~isempty(st.vols{volhandle}.roi.xyz)
      st.vols{volhandle}.roi.xyz = union(st.vols{volhandle}.roi.xyz',toset','rows')';
    else
      st.vols{volhandle}.roi.xyz = toset;
    end;
  end;

  if isfield(st.vols{volhandle}, 'blobs')
    nblobs=length(st.vols{volhandle}.blobs);
    if nblobs>1
      blobstmp(1:st.vols{volhandle}.roi.hroi-1) = st.vols{volhandle}.blobs(1:st.vols{volhandle}.roi.hroi-1);
      blobstmp(st.vols{volhandle}.roi.hroi:nblobs-1) = st.vols{volhandle}.blobs(st.vols{volhandle}.roi.hroi+1:nblobs);
      st.vols{volhandle}.blobs=blobstmp;
    else
      rmfield(st.vols{volhandle},'blobs');
    end;
  end;
  if isfield(st.vols{volhandle}, 'blobs')
    st.vols{volhandle}.roi.hroi = prod(size(st.vols{volhandle}.blobs))+1;
  else
    st.vols{volhandle}.roi.hroi = 1;
  end;
  if isempty(st.vols{volhandle}.roi.xyz)  % initialised with empty roi
    spm_orthviews('addcolouredblobs', volhandle, ...
      [1; 1; 1], 0, st.vols{volhandle}.roi.Vroi.mat,[1 3 1]); 
  else
    spm_orthviews('addcolouredblobs', volhandle, ...
      st.vols{volhandle}.roi.xyz, ones(size(st.vols{volhandle}.roi.xyz,2),1), ... 
      st.vols{volhandle}.roi.Vroi.mat,[1 3 1]); % use color that is more intense than standard rgb range
  end;
  st.vols{volhandle}.blobs{st.vols{volhandle}.roi.hroi}.max=2;
  spm_orthviews('redraw');
end;

spm('pointer','arrow');


function varargout = stack(cmd, varargin)

switch cmd
  case 'init'
    varargout{1}.st = cell(varargin{1},1);
    varargout{1}.top= 0;

  case 'isempty'
    varargout{1} = (varargin{1}.top==0);

  case 'push'
    stck = varargin{1};
    if (stck.top < size(stck.st,1))
      stck.top = stck.top + 1;
      stck.st{stck.top} = varargin{2};
      varargout{1}=stck;
    else
      error('Stack overflow\n');
    end;

  case 'pop'
    if stack('isempty',varargin{1})
      error('Stack underflow\n');
    else
      varargout{2} = varargin{1}.st{varargin{1}.top};
      varargin{1}.top = varargin{1}.top - 1;
      varargout{1} = varargin{1};
    end;
end;

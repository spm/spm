function ret = spm_ov_quiver3d(varargin)
% function ret = spm_ov_quiver3d(varargin)
%
% NAME   quiver3d    (Plugin to spm_orthviews)
%
% FORMAT ret = spm_orthviews('quiver3d', command, volhandle[, varargin])
%
% COMMANDS
% spm_orthviews('quiver3d', 'init', volhandle, Ve1fnames, Ve2fnames,
%                                   Ve3fnames, Vlafnames, Vmaskfname 
%                                   [, maskmin, thresh(2), qst, ql])
%     ARGUMENTS
%       Ve1fnames   Eigenvector images (3 per eigenvector e1, e2, e3)
%       Ve2fnames
%       Ve3fnames
%       Vlafnames   Eigenvalue images (1 per eigenvector)
%       Vmaskfname  Mask image (only quivers within mask will be displayed)
%     [optional]
%       thresh      (Default [.1 Inf]) Thresholds for masking
%       qst         (Default 3)  Stepsize for quiver
%       ql          (Default 5)  Size for quiver
%
% spm_orthviews('quiver3d', 'redraw', volhandle)
%
% spm_orthviews('quiver3d, 'delete', volhandle)
%_____________________________________________________________________________
% %W% Volkmar Glauche <glauche@uke.uni-hamburg.de> %E%

global st;
if isempty(st)
  error('spm_ov_quiver3d: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
  error('spm_orthviews(''quiver3d'',cmd,volhandle): Not enough arguments');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};
switch cmd
  case 'init'
    %function addquiver3d(handle, Ve1fnames, Ve2fnames, Ve3fnames, Vlafnames, Vmaskfname, varargin)
    if nargin < 7
      error('spm_orthviews(''quiver3d'', ''init'',...): Not enough arguments');
    end;

    Ve1 = spm_vol(varargin{3});
    Ve2 = spm_vol(varargin{4});
    Ve3 = spm_vol(varargin{5});
    Vla = spm_vol(varargin{6});
    Vmask = spm_vol(varargin{7});
    if all([length(Ve1) length(Ve2) length(Ve3) length(Vla)] == 3)
      st.vols{volhandle}.quiver3d = struct('e1',Ve1,'e2',Ve2,'e3',Ve3, 'la', Vla, ...
	  'mask',Vmask, 'thresh', [.1 Inf], ...
	  'qst',3,'ql', 5,'qht',[],'qhc',[],'qhs',[]);
    else
      error('spm_orthviews(''quiver3d'', ''init'',...): Please specify 3 images for each eigenvector!');
    end;
    if nargin > 7
      if ~isempty(varargin{8})
	st.vols{volhandle}.quiver.thresh(1) = varargin{8};
      end;
    end;
    if nargin > 8
      if ~isempty(varargin{9})
	st.vols{volhandle}.quiver.thresh(2) = varargin{9};
      end;
    end
    if nargin > 9
      if ~isempty(varargin{10})
	st.vols{volhandle}.quiver.qst = varargin{10};
      end;
    end;
    if nargin > 10
      if ~isempty(varargin{11})
	st.vols{volhandle}.quiver.ql = varargin{11};
      end;
    end;
    
  case 'redraw'
    TM0 = varargin{3};
    TD  = varargin{4};
    CM0 = varargin{5};
    CD  = varargin{6};
    SM0 = varargin{7};
    SD  = varargin{8};
    if isfield(st.vols{volhandle},'quiver3d')
      q3d = st.vols{volhandle}.quiver3d;
      % need to delete old quiver lines before redrawing
      if ishandle(q3d.qht)
	delete(q3d.qht);
      end;
      if ishandle(q3d.qhc)
	delete(q3d.qhc);
      end;
      if ishandle(q3d.qhs)
	delete(q3d.qhs);
      end;
      
      % step size for selection of locations
      prm = spm_imatrix(st.Space);
      qst = ceil(q3d.qst/prm(7));
      qst1 = ceil(qst/2);
      
      Me1x   = st.Space\(st.vols{volhandle}.premul*q3d.e1(1).mat);
      Me1y   = st.Space\(st.vols{volhandle}.premul*q3d.e1(2).mat);
      Me1z   = st.Space\(st.vols{volhandle}.premul*q3d.e1(3).mat);
      Me2x   = st.Space\(st.vols{volhandle}.premul*q3d.e2(1).mat);
      Me2y   = st.Space\(st.vols{volhandle}.premul*q3d.e2(2).mat);
      Me2z   = st.Space\(st.vols{volhandle}.premul*q3d.e2(3).mat);
      Me3x   = st.Space\(st.vols{volhandle}.premul*q3d.e3(1).mat);
      Me3y   = st.Space\(st.vols{volhandle}.premul*q3d.e3(2).mat);
      Me3z   = st.Space\(st.vols{volhandle}.premul*q3d.e3(3).mat);
      Mla1   = st.Space\(st.vols{volhandle}.premul*q3d.la(1).mat);
      Mla2   = st.Space\(st.vols{volhandle}.premul*q3d.la(2).mat);
      Mla3   = st.Space\(st.vols{volhandle}.premul*q3d.la(2).mat);
      Mm     = st.Space\(st.vols{volhandle}.premul*q3d.mask.mat);

      % transversal: Ft, Vt, FVCt
      % laX, eX, x, y, xyz, mask will be reused for all 3 views
      la1	 = spm_slice_vol(q3d.la(1),inv(TM0*Mla1),TD,0);
      la1(la1==0) = NaN;
      la2	 = spm_slice_vol(q3d.la(2),inv(TM0*Mla2),TD,0)./la1; % scale to la1
      la3	 = spm_slice_vol(q3d.la(3),inv(TM0*Mla3),TD,0)./la1; % scale to la1
      e1   = cat(3, spm_slice_vol(q3d.e1(1),inv(TM0*Me1x),TD,0), ...
	  spm_slice_vol(q3d.e1(2),inv(TM0*Me1y),TD,0), ...
	  spm_slice_vol(q3d.e1(3),inv(TM0*Me1z),TD,0));
      e2   = cat(3, la2.*spm_slice_vol(q3d.e2(1),inv(TM0*Me2x),TD,0), ...
	  la2.*spm_slice_vol(q3d.e2(2),inv(TM0*Me2y),TD,0), ...
	  la2.*spm_slice_vol(q3d.e2(3),inv(TM0*Me2z),TD,0));
      e3   = cat(3, la3.*spm_slice_vol(q3d.e3(1),inv(TM0*Me3x),TD,0), ...
	  la3.*spm_slice_vol(q3d.e3(2),inv(TM0*Me3y),TD,0), ...
	  la3.*spm_slice_vol(q3d.e3(3),inv(TM0*Me3z),TD,0));
      mask = spm_slice_vol(q3d.mask,inv(TM0*Mm),TD,0);
      mask = ((mask > q3d.thresh(1)) & (mask < q3d.thresh(2)))|isnan(la1);
      e1   = e1(qst1:qst:end,qst1:qst:end,:);
      e2   = e2(qst1:qst:end,qst1:qst:end,:);
      e3   = e3(qst1:qst:end,qst1:qst:end,:);
      mask = mask(qst1:qst:end,qst1:qst:end);
      nmask = sum(mask(:));
      se = size(e1);

      % rotate vectors according to the plane's coordinate system			 
      e1 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e1,prod(se(1:2)),3)';
      e2 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e2,prod(se(1:2)),3)';
      e3 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e3,prod(se(1:2)),3)';
      
      [x y] = ndgrid([1:TD(1)]-.5,[1:TD(2)]-.5);
      x = x(qst1:qst:end,qst1:qst:end);
      y = y(qst1:qst:end,qst1:qst:end);
      xyz   = [x(mask(:)), y(mask(:)), zeros(nmask,1)];

      % create vertices
      Vt = cat(1,xyz+q3d.ql*e1(:,mask(:))',...
	  xyz-q3d.ql*e1(:,mask(:))',...
	  xyz+q3d.ql*e2(:,mask(:))',...
	  xyz-q3d.ql*e2(:,mask(:))',...
	  xyz+q3d.ql*e3(:,mask(:))',...
	  xyz-q3d.ql*e3(:,mask(:))');

      % create faces
      Ft = [repmat([1:nmask],1,4) nmask+repmat([1:nmask],1,4);...
	repmat([repmat(2*nmask+[1:nmask],1,2) ...
	      repmat(3*nmask+[1:nmask],1,2)],1,2);...
	repmat([4*nmask+[1:nmask] 5*nmask+[1:nmask]],1,4)]';

      % create colours
      FVCt = cat(1,repmat([1 0 0],2*nmask,1),...
	  repmat([0 1 0],2*nmask,1),...
	  repmat([0 0 1],2*nmask,1));

      % coronal: Fc, Vc, FVCc
      % laX, eX, x, y, xyz will be reused for all 3 views
      la1	 = spm_slice_vol(q3d.la(1),inv(CM0*Mla1),CD,0);
      la1(la1==0) = NaN;
      la2	 = spm_slice_vol(q3d.la(2),inv(CM0*Mla2),CD,0)./la1; % scale to la1
      la3	 = spm_slice_vol(q3d.la(3),inv(CM0*Mla3),CD,0)./la1; % scale to la1
      e1   = cat(3, spm_slice_vol(q3d.e1(1),inv(CM0*Me1x),CD,0), ...
	  spm_slice_vol(q3d.e1(2),inv(CM0*Me1y),CD,0), ...
	  spm_slice_vol(q3d.e1(3),inv(CM0*Me1z),CD,0));
      e2   = cat(3, la2.*spm_slice_vol(q3d.e2(1),inv(CM0*Me2x),CD,0), ...
	  la2.*spm_slice_vol(q3d.e2(2),inv(CM0*Me2y),CD,0), ...
	  la2.*spm_slice_vol(q3d.e2(3),inv(CM0*Me2z),CD,0));
      e3   = cat(3, la3.*spm_slice_vol(q3d.e3(1),inv(CM0*Me3x),CD,0), ...
	  la3.*spm_slice_vol(q3d.e3(2),inv(CM0*Me3y),CD,0), ...
	  la3.*spm_slice_vol(q3d.e3(3),inv(CM0*Me3z),CD,0));
      mask = spm_slice_vol(q3d.mask,inv(CM0*Mm),CD,0);
      mask = ((mask > q3d.thresh(1)) & (mask < q3d.thresh(2)))|isnan(la1);
      e1   = e1(qst1:qst:end,qst1:qst:end,:);
      e2   = e2(qst1:qst:end,qst1:qst:end,:);
      e3   = e3(qst1:qst:end,qst1:qst:end,:);
      mask = mask(qst1:qst:end,qst1:qst:end);
      nmask = sum(mask(:));
      se = size(e1);

      % rotate vectors according to the plane's coordinate system			 
      e1 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e1,prod(se(1:2)),3)';
      e2 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e2,prod(se(1:2)),3)';
      e3 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e3,prod(se(1:2)),3)';
      
      % sort eigenvector components (coronal - x z y)
      e1 = e1([1 3 2],:);
      e2 = e2([1 3 2],:);
      e3 = e3([1 3 2],:);
      
      [x y] = ndgrid([1:CD(1)]-.5,[1:CD(2)]-.5);
      x = x(qst1:qst:end,qst1:qst:end);
      y = y(qst1:qst:end,qst1:qst:end);
      xyz   = [x(mask(:)), y(mask(:)), zeros(nmask,1)];

      % create vertices
      Vc = cat(1,xyz+q3d.ql*e1(:,mask(:))',...
	  xyz-q3d.ql*e1(:,mask(:))',...
	  xyz+q3d.ql*e2(:,mask(:))',...
	  xyz-q3d.ql*e2(:,mask(:))',...
	  xyz+q3d.ql*e3(:,mask(:))',...
	  xyz-q3d.ql*e3(:,mask(:))');

      % create faces
      Fc = [repmat([1:nmask],1,4) nmask+repmat([1:nmask],1,4);...
	repmat([repmat(2*nmask+[1:nmask],1,2) ...
	      repmat(3*nmask+[1:nmask],1,2)],1,2);...
	repmat([4*nmask+[1:nmask] 5*nmask+[1:nmask]],1,4)]';

      % create colours
      FVCc = cat(1,repmat([1 0 0],2*nmask,1),...
	  repmat([0 1 0],2*nmask,1),...
	  repmat([0 0 1],2*nmask,1));
      
      % sagittal: Fs, Vs, FVCs
      % laX, eX, x, y, xyz will be reused for all 3 views
      la1	 = spm_slice_vol(q3d.la(1),inv(SM0*Mla1),SD,0);
      la1(la1==0) = NaN;
      la2	 = spm_slice_vol(q3d.la(2),inv(SM0*Mla2),SD,0)./la1; % scale to la1
      la3	 = spm_slice_vol(q3d.la(3),inv(SM0*Mla3),SD,0)./la1; % scale to la1
      e1   = cat(3, spm_slice_vol(q3d.e1(1),inv(SM0*Me1x),SD,0), ...
	  spm_slice_vol(q3d.e1(2),inv(SM0*Me1y),SD,0), ...
	  spm_slice_vol(q3d.e1(3),inv(SM0*Me1z),SD,0));
      e2   = cat(3, la2.*spm_slice_vol(q3d.e2(1),inv(SM0*Me2x),SD,0), ...
	  la2.*spm_slice_vol(q3d.e2(2),inv(SM0*Me2y),SD,0), ...
	  la2.*spm_slice_vol(q3d.e2(3),inv(SM0*Me2z),SD,0));
      e3   = cat(3, la3.*spm_slice_vol(q3d.e3(1),inv(SM0*Me3x),SD,0), ...
	  la3.*spm_slice_vol(q3d.e3(2),inv(SM0*Me3y),SD,0), ...
	  la3.*spm_slice_vol(q3d.e3(3),inv(SM0*Me3z),SD,0));
      mask = spm_slice_vol(q3d.mask,inv(SM0*Mm),SD,0);
      mask = ((mask > q3d.thresh(1)) & (mask < q3d.thresh(2)))|isnan(la1);
      e1   = e1(qst1:qst:end,qst1:qst:end,:);
      e2   = e2(qst1:qst:end,qst1:qst:end,:);
      e3   = e3(qst1:qst:end,qst1:qst:end,:);
      mask = mask(qst1:qst:end,qst1:qst:end);
      nmask = sum(mask(:));
      se = size(e1);

      % rotate vectors according to the plane's coordinate system			 
      e1 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e1,prod(se(1:2)),3)';
      e2 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e2,prod(se(1:2)),3)';
      e3 = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(e3,prod(se(1:2)),3)';
      
      % sort eigenvector components (sagittal - -y z x)
      e1 = [-e1(2,:); e1([3 1],:)];
      e2 = [-e2(2,:); e2([3 1],:)];
      e3 = [-e2(2,:); e3([3 1],:)];
      
      [x y] = ndgrid([1:SD(1)]-.5,[1:SD(2)]-.5);
      x = x(qst1:qst:end,qst1:qst:end);
      y = y(qst1:qst:end,qst1:qst:end);
      xyz   = [x(mask(:)), y(mask(:)), zeros(nmask,1)];

      % create vertices
      Vs = cat(1,xyz+q3d.ql*e1(:,mask(:))',...
	  xyz-q3d.ql*e1(:,mask(:))',...
	  xyz+q3d.ql*e2(:,mask(:))',...
	  xyz-q3d.ql*e2(:,mask(:))',...
	  xyz+q3d.ql*e3(:,mask(:))',...
	  xyz-q3d.ql*e3(:,mask(:))');

      % create faces
      Fs = [repmat([1:nmask],1,4) nmask+repmat([1:nmask],1,4);...
	repmat([repmat(2*nmask+[1:nmask],1,2) ...
	      repmat(3*nmask+[1:nmask],1,2)],1,2);...
	repmat([4*nmask+[1:nmask] 5*nmask+[1:nmask]],1,4)]';

      % create colours
      FVCs = cat(1,repmat([1 0 0],2*nmask,1),...
	  repmat([0 1 0],2*nmask,1),...
	  repmat([0 0 1],2*nmask,1));
      
      % transversal - plot (x y z)
      np = get(st.vols{volhandle}.ax{1}.ax,'NextPlot');
      set(st.vols{volhandle}.ax{1}.ax,'NextPlot','add');
      axes(st.vols{volhandle}.ax{1}.ax);
      st.vols{volhandle}.quiver3d.qht = patch('Faces', Ft, 'Vertices',Vt,...
	  'FaceVertexCData', FVCt, ...
	  'Edgecolor', 'interp', 'FaceColor', 'interp');
      set(st.vols{volhandle}.ax{1}.ax,'NextPlot',np);
      set(st.vols{volhandle}.quiver3d.qht, 'Parent',st.vols{volhandle}.ax{1}.ax, 'HitTest','off');
      
      % coronal - plot (x z y)
      np = get(st.vols{volhandle}.ax{2}.ax,'NextPlot');
      set(st.vols{volhandle}.ax{2}.ax,'NextPlot','add');
      axes(st.vols{volhandle}.ax{2}.ax);
      st.vols{volhandle}.quiver3d.qhc = patch('Faces', Fc, 'Vertices',Vc,...
	  'FaceVertexCData', FVCc, ...
	  'Edgecolor', 'interp', 'FaceColor', 'interp'); 
      set(st.vols{volhandle}.ax{2}.ax,'NextPlot',np);
      set(st.vols{volhandle}.quiver3d.qhc, 'Parent',st.vols{volhandle}.ax{2}.ax, 'HitTest','off');
      
      % sagittal - plot (-y z x)
      np = get(st.vols{volhandle}.ax{3}.ax,'NextPlot');
      set(st.vols{volhandle}.ax{3}.ax,'NextPlot','add');
      axes(st.vols{volhandle}.ax{3}.ax);
      st.vols{volhandle}.quiver3d.qhs = patch('Faces', Fs, 'Vertices',Vs,...
	  'FaceVertexCData', FVCs, ...
	  'Edgecolor', 'interp', 'FaceColor', 'interp');
      set(st.vols{volhandle}.ax{3}.ax,'NextPlot',np);
      set(st.vols{volhandle}.quiver3d.qhs, 'Parent',st.vols{volhandle}.ax{3}.ax, 'HitTest','off');
    end; %quiver3d

  case 'delete'
    if isfield(st.vols{volhandle},'quiver3d'),
      delete(st.vols{volhandle}.quiver3d.qht);
      delete(st.vols{volhandle}.quiver3d.qhc);
      delete(st.vols{volhandle}.quiver3d.qhs);
      st.vols{volhandle} = rmfield(st.vols{volhandle},'quiver3d');
    end;
  %-------------------------------------------------------------------------
  % Context menu and callbacks
  case 'context_menu'  
    item0 = uimenu(varargin{3}, 'Label', 'Quiver3D');
      item1 = uimenu(item0, 'Label', 'Add', 'Callback', ...
	  ['feval(''spm_ov_quiver3d'',''context_init'', ', ...
	    num2str(volhandle), ');'], 'Tag', ['QUIVER3D_0_', num2str(volhandle)]);
      item2 = uimenu(item0, 'Label', 'Properties', ...
	  'Visible', 'off', 'Tag', ['QUIVER3D_1_', num2str(volhandle)]);
        item2_1 = uimenu(item2, 'Label', 'Mask threshold', 'Callback', ...
	    ['feval(''spm_ov_quiver3d'',''context_edit'',', ...
	      num2str(volhandle), ',''thresh'');']);
	item2_2 = uimenu(item2, 'Label', 'Quiver3D distance', 'Callback', ...
	    ['feval(''spm_ov_quiver3d'',''context_edit'',', ...
	      num2str(volhandle), ',''qst'');']); 
	item2_3 = uimenu(item2, 'Label', 'Quiver3D length', 'Callback', ...
	    ['feval(''spm_ov_quiver3d'',''context_edit'',', ...
	      num2str(volhandle), ',''ql'');']); 
      item3 = uimenu(item0, 'Label', 'Remove', 'Callback', ...
	  ['feval(''spm_ov_quiver3d'',''context_delete'', ', ...
	    num2str(volhandle), ');'], 'Visible', 'off', ...
	  'Tag', ['QUIVER3D_1_', num2str(volhandle)]);

  case 'context_init'
    Finter = spm_figure('FindWin', 'Interactive');
    spm_input('!DeleteInputObj',Finter);
    Ve1fnames = spm_get(3,'evec1*.img','Components of 1st eigenvector');
    Ve2fnames = spm_get(3,'evec2*.img','Components of 2nd eigenvector');
    Ve3fnames = spm_get(3,'evec3*.img','Components of 3rd eigenvector');
    Vlafnames = spm_get(3,'eval*.img','Eigenvalue images');
    Vmaskfname = spm_get(1,'*.img','Mask image');
    feval('spm_ov_quiver3d','init',volhandle,Ve1fnames,Ve2fnames,Ve3fnames,Vlafnames,Vmaskfname);
    obj = findobj(0, 'Tag',  ['QUIVER3D_1_', num2str(volhandle)]);
    set(obj, 'Visible', 'on');
    obj = findobj(0, 'Tag',  ['QUIVER3D_0_', num2str(volhandle)]);
    set(obj, 'Visible', 'off');
    spm_orthviews('redraw');
  
  case 'context_edit'
    Finter = spm_figure('FindWin', 'Interactive');
    spm_input('!DeleteInputObj',Finter);
    switch varargin{3}
      case 'thresh'
	in = spm_input('Mask threshold {min max}','!+1','e', ...
	    num2str(st.vols{volhandle}.quiver3d.thresh), [1 2]);
      case 'qst'
	in = spm_input('Quiver3D distance','!+1','e', ...
	    num2str(st.vols{volhandle}.quiver3d.qst), 1);
      case 'ql'
	in = spm_input('Quiver3D length','!+1','e', ...
	    num2str(st.vols{volhandle}.quiver3d.ql), 1);
    end;
    spm_input('!DeleteInputObj',Finter);
    st.vols{volhandle}.quiver3d = setfield(st.vols{volhandle}.quiver3d, ...
	varargin{3}, in);
    spm_orthviews('redraw');
    
  case 'context_delete'
    feval('spm_ov_quiver3d','delete',volhandle);
    obj = findobj(0, 'Tag',  ['QUIVER3D_1_', num2str(volhandle)]);
    set(obj, 'Visible', 'off');
    obj = findobj(0, 'Tag',  ['QUIVER3D_0_', num2str(volhandle)]);
    set(obj, 'Visible', 'on');
    
  otherwise
    warning('spm_orthviews(''quiver3d'',...): Unknown action string');
  end;


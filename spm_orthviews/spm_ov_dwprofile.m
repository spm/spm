function ret = spm_ov_dwprofile(varargin)
% This routine is a plugin to spm_orthviews for SPM5. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the matlab prompt.
%_______________________________________________________________________
%
% @(#) $Id: spm_ov_dwprofile.m,v 1.8 2006/02/21 13:58:09 glauche Exp $

global st;
if isempty(st)
  error('dwprofile: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
  error(['dwprofile: Wrong number of arguments. Usage: spm_orthviews(''dwprofile'', ' ...
	 'cmd, volhandle, varargin)']);
end;

cmd = lower(varargin{1});
volhandle = varargin{2};
switch cmd
 case 'init'
  if isstruct(varargin{3})
    st.vols{volhandle}.dwprofile=varargin{3};
  else
    voi = varargin{3};
    dwprofile=dti_surface('defaults');
    dwprofile.scaling=0;
    dwprofile.res=.2;
    for k=1:numel(voi)
      dwprofile.voi=lower(voi{k});
      st.vols{volhandle}.dwprofile(k)=dwprofile;
    end;
  end;
  spm_orthviews('redraw');
  for k = 1:numel(st.vols{volhandle}.dwprofile)
    if strcmp(st.vols{volhandle}.dwprofile(k).voi(1:5),'plane')
      spm_ovhelper_3Dreg('register', ...
			 st.vols{volhandle}.dwprofile(k).f,...
			 st.vols{volhandle});
    end;
    set(st.vols{volhandle}.dwprofile(k).f, 'CloseRequestFcn',...
		      ['feval(''' mfilename ''', ''context_delete'',' ...
		       num2str(volhandle) ');']);
  end;
  return;
 case 'redraw'
  TM0 = varargin{3};
  TD  = varargin{4};
  CM0 = varargin{5};
  CD  = varargin{6};
  SM0 = varargin{7};
  SD  = varargin{8};

  for k=1:numel(st.vols{volhandle}.dwprofile)
    st.vols{volhandle}.dwprofile(k).xyzvx=[];
    voi{k} = lower(st.vols{volhandle}.dwprofile(k).voi);
    pos = spm_orthviews('pos')';
    vx = sqrt(sum(st.Space(1:3,1:3).^2));
    switch voi{k}
     case 'point'
      st.vols{volhandle}.dwprofile(k).xyzmm=pos;
     case 'planexy'
      [x y z]=ndgrid(vx(1)*[st.bb(1,1):st.bb(2,1)], ...
                     vx(2)*[st.bb(1,2):st.bb(2,2)], pos(3));
      st.vols{volhandle}.dwprofile(k).xyzmm = [x(:), y(:), z(:)];
     case 'planexz'
      [x y z]=ndgrid(vx(1)*[st.bb(1,1):st.bb(2,1)], pos(2), ...
                     vx(3)*[st.bb(1,3):st.bb(2,3)]);
      st.vols{volhandle}.dwprofile(k).xyzmm = [x(:), y(:), z(:)];
     case 'planeyz'
      [x y z]=ndgrid(pos(1), vx(2)*[st.bb(1,2):st.bb(2,2)], ...
                         vx(3)*[st.bb(1,3):st.bb(2,3)]);
      st.vols{volhandle}.dwprofile(k).xyzmm = [x(:), y(:), z(:)];
    end;
  end;
  [res st.vols{volhandle}.dwprofile]= ...
      dti_surface(st.vols{volhandle}.dwprofile);
  for k=1:numel(st.vols{volhandle}.dwprofile)
    st.vols{volhandle}.dwprofile(k).voi = voi{k};% keep voi setting
    if strncmp(st.vols{volhandle}.dwprofile(k).voi,'plane',5) 
      set(st.vols{volhandle}.dwprofile(k).ax, 'ButtonDownFcn', ...
                        sprintf('spm_ov_dwprofile(''setpos'',''%s'')', ...
                                st.vols{volhandle}.dwprofile(k).voi(6:7)));
    end;
  end;
 
 case 'delete'
  if isfield(st.vols{volhandle},'dwprofile'),
    delete([st.vols{volhandle}.dwprofile(:).f]);
    st.vols{volhandle} = rmfield(st.vols{volhandle},'dwprofile');
  end;
  %-------------------------------------------------------------------------
  % Context menu and callbacks
 case 'context_menu'  
  item0 = uimenu(varargin{3}, 'Label', 'DW Profile');
  item1a = uimenu(item0, 'Label', 'Current voxel', 'Callback', ...
                 ['feval(''spm_ov_dwprofile'',''context_init'', ', ...
                  num2str(volhandle), ',''point'');'], ...
                 'Tag',['DWPROFILE_0_', num2str(volhandle)]);
  item1b = uimenu(item0, 'Label', 'XYZ planes', 'Callback', ...
                 ['feval(''spm_ov_dwprofile'',''context_init'', ', ...
                  num2str(volhandle), ',''planes'');'], ...
                 'Tag',['DWPROFILE_0_', num2str(volhandle)]);
  item1c = uimenu(item0, 'Label', 'Both', 'Callback', ...
                 ['feval(''spm_ov_dwprofile'',''context_init'', ', ...
                  num2str(volhandle), ',''both'');'], ...
                 'Tag',['DWPROFILE_0_', num2str(volhandle)]);
  item2 = uimenu(item0, 'Label', 'Properties', ...
                 'Visible', 'off', 'Tag', ['DWPROFILE_1_', num2str(volhandle)]);
  item2_1 = uimenu(item2, 'Label', 'Sphere resolution', 'Callback', ...
                   ['feval(''spm_ov_dwprofile'',''context_edit'',', ...
                    num2str(volhandle), ',''res'');']);
  item2_2 = uimenu(item2, 'Label', 'Transparency');
  item2_2_1 = uimenu(item2_2, 'Label','On', 'Callback', ...
		     ['feval(''spm_ov_dwprofile'',''context_edit'',', ...
		      num2str(volhandle), ',''alphaon'');'],...
		     'Tag', ['DWPROFILE_1_ALPHAON', num2str(volhandle)]);
  item2_2_2 = uimenu(item2_2, 'Label', 'Off', 'Checked','On', 'Callback', ...
		     ['feval(''spm_ov_dwprofile'',''context_edit'',', ...
		      num2str(volhandle), ',''alphaoff'');'],...
		     'Tag', ['DWPROFILE_1_ALPHAOFF', num2str(volhandle)]);
  % 	item2_3 = uimenu(item2, 'Label', 'Dwprofile distance', 'Callback', ...
  % 	    ['feval(''spm_ov_dwprofile'',''context_edit'',', ...
  % 	      num2str(volhandle), ',''qst'');']); 
  % 	item2_4 = uimenu(item2, 'Label', 'Dwprofile length', 'Callback', ...
  % 	    ['feval(''spm_ov_dwprofile'',''context_edit'',', ...
  % 	      num2str(volhandle), ',''ql'');']); 
  % 	item2_5 = uimenu(item2, 'Label', 'Linewidth', 'Callback', ...
  % 	    ['feval(''spm_ov_dwprofile'',''context_edit'',', ...
  % 	      num2str(volhandle), ',''qw'');']); 
  item3 = uimenu(item0, 'Label', 'Remove', 'Callback', ...
                 ['feval(''spm_ov_dwprofile'',''context_delete'', ', ...
                  num2str(volhandle), ');'], 'Visible', 'off', ...
                 'Tag', ['DWPROFILE_1_', num2str(volhandle)]);
  
 case 'context_init'
  Finter = spm_figure('FindWin', 'Interactive');
  spm_input('!DeleteInputObj',Finter);
  dwprofile=dti_surface('defaults');
  dwprofile.scaling=0;
  dwprofile.alpha=0;
  nd=spm_input('# of datasets','!+1','n',1,1);
  for d=1:nd
    dwprofile.order = spm_input(sprintf('Set %d - Tensor order',d), 1,'n', ...
                               2,1);
    nimg = (dwprofile.order+1)*(dwprofile.order+2)/2;
    selstr = sprintf('D%s_.*',repmat('[xyz]',1,dwprofile.order));
    dwprofile.files = cellstr(spm_select(nimg,'image',...
                              sprintf('Set %d -DT images',d),'',pwd,selstr));
    switch lower(varargin{3})
     case 'point'
      voi={'point'};
     case 'planes'
      voi={'planexy','planexz','planeyz'};
     case 'both'
      voi={'planexy','planexz','planeyz','point'};
    end;      
    for k=1:numel(voi)
      dwprofile.voi=lower(voi{k});
      if strcmp(dwprofile.voi,'point')
        dwprofile.res=.05;
      else
        dwprofile.res=.2;
      end;
      fdwprofile((d-1)*numel(voi)+k)=dwprofile;
    end;
  end;  
  spm_input('!DeleteInputObj',Finter);
  feval('spm_ov_dwprofile','init',volhandle,fdwprofile);
  obj = findobj(0, 'Tag',  ['DWPROFILE_1_', num2str(volhandle)]);
  set(obj, 'Visible', 'on');
  obj = findobj(0, 'Tag',  ['DWPROFILE_0_', num2str(volhandle)]);
  set(obj, 'Visible', 'off');
  
 case 'context_edit'
  Finter = spm_figure('FindWin', 'Interactive');
  spm_input('!DeleteInputObj',Finter);
  switch varargin{3}
   case 'res'
    points=find(strncmp(cellstr(strvcat(st.vols{volhandle}.dwprofile.voi)), ...
                   'point',5));
    planes=find(strncmp(cellstr(strvcat(st.vols{volhandle}.dwprofile.voi)), ...
                   'plane',5));
    if ~isempty(points)
      in = spm_input('Single points: sphere resolution','!+1','e', ...
                     num2str(st.vols{volhandle}.dwprofile(points(1)).res), ...
                     1);
      [st.vols{volhandle}.dwprofile(points).(varargin{3})] = deal(in);
    end;
    if ~isempty(planes)
      in = spm_input('Planes: sphere resolution','!+1','e', ...
                     num2str(st.vols{volhandle}.dwprofile(planes(1)).res), ...
                     1);
      [st.vols{volhandle}.dwprofile(planes).(varargin{3})] = deal(in);
    end;
   case {'alphaon','alphaoff'}
    halphaon = findobj(get(gcbo,'Parent'), 'Tag',['DWPROFILE_1_ALPHAON', ...
		    num2str(volhandle)]);
    halphaoff = findobj(get(gcbo,'Parent'), 'Tag',['DWPROFILE_1_ALPHAOFF', ...
		    num2str(volhandle)]);
    switch(varargin{3}(6:end))
     case 'on'
      set(halphaon,'Checked','on');
      set(halphaoff,'Checked','off');
      [st.vols{volhandle}.dwprofile.alpha] = deal(1);
     case 'off'
      set(halphaon,'Checked','off');
      set(halphaoff,'Checked','on');
      [st.vols{volhandle}.dwprofile.alpha] = deal(0);
    end;
  end;
  spm_input('!DeleteInputObj',Finter);
  spm_orthviews('redraw');
    
 case 'context_delete'
  % pass through optional 3rd argument
  feval('spm_ov_dwprofile','delete',volhandle,varargin{3:end});
  obj = findobj(0, 'Tag',  ['DWPROFILE_1_', num2str(volhandle)]);
  set(obj, 'Visible', 'off');
  obj = findobj(0, 'Tag',  ['DWPROFILE_0_', num2str(volhandle)]);
  set(obj, 'Visible', 'on');

 case 'setpos'
  axpos=(get(gcbo,'CurrentPoint'));
  ovpos=spm_orthviews('pos');
  switch varargin{2}
   case 'xy'
    t = (ovpos(3)-axpos(2,3))/(axpos(1,3)-axpos(2,3));
   case 'xz'
    t = (ovpos(2)-axpos(2,2))/(axpos(1,2)-axpos(2,2));
   case 'yz'
    t = (ovpos(1)-axpos(2,1))/(axpos(1,1)-axpos(2,1));
  end;
  if t>=0 & t<=1
    spm_orthviews('reposition',[t*(axpos(1,1)-axpos(2,1))+axpos(2,1), ...
                        t*(axpos(1,2)-axpos(2,2))+axpos(2,2),...
                        t*(axpos(1,3)-axpos(2,3))+axpos(2,3)]);
  end;
 otherwise

  fprintf('spm_orthviews(''dwprofile'',...): Unknown action string %s', cmd);
end;

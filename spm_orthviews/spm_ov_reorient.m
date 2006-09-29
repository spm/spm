function ret = spm_ov_reorient(varargin)
% Reorient tool - plugin for spm_orthviews
%
% This tool provides the capabilities of the reorientation widget in SPMs
% "DISPLAY" for any image displayed within spm_orthviews. The control fields
% are drawn in the SPM interactive window and work as described in the
% Display routine. 
% The advantage of using this tool within CheckReg is that it allows to
% reorient images while comparing their position to reference images
% simultaneously. 
%
% This routine is a plugin to spm_orthviews for SPM5. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the matlab prompt.
%_____________________________________________________________________________
% $Id: spm_ov_reorient.m 640 2006-09-29 09:53:58Z volkmar $

rev = '$Revision: 640 $';

global st;
if isempty(st)
  error('reorient: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
  error('reorient: Wrong number of arguments. Usage: spm_orthviews(''reorient'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};
switch cmd
  %-------------------------------------------------------------------------
  % Context menu and callbacks
 case 'context_menu'  
  item0 = uimenu(varargin{3}, 'Label', 'Reorient image(s)',...
		 'Tag', ['REORIENT_M_', num2str(volhandle)]);
  item1 = uimenu(item0, 'Label', 'All images', 'Callback', ...
		 ['spm_orthviews(''reorient'',''context_init'', 0);'],...
		 'Tag', ['REORIENT_0_', num2str(volhandle)]);
  item2 = uimenu(item0, 'Label', 'Current image', 'Callback', ...
		 ['spm_orthviews(''reorient'',''context_init'', ', ...
		  num2str(volhandle), ');'],...
		 'Tag', ['REORIENT_0_', num2str(volhandle)]);
  item3 = uimenu(item0, 'Label', 'Quit Reorient image', ...
		 'Tag', ['REORIENT_1_', num2str(volhandle)], ...
		 'Visible', 'off');
  item1 = uimenu(item0, 'Label', 'Help', 'Callback', ...
                 ['feval(''spm_help'',''' mfilename ''');']);
  ret = item0;
  
 case 'context_init'
  Finter = spm_figure('FindWin', 'Interactive');
  Fgraph = spm_figure('FindWin', 'Graphics');
  figure(Finter);
  spm_input('!DeleteInputObj',Finter);
  handles = spm_orthviews('valid_handles');
  labels = {'right  {mm}', 'forward  {mm}', 'up  {mm}',...
	    'pitch  {rad}', 'roll  {rad}', 'yaw  {rad}',...
	    'resize  {x}', 'resize  {y}', 'resize {z}'};
  tooltips = {'translate', 'translate', 'translate', 'rotate', 'rotate', ...
	      'rotate', 'zoom', 'zoom', 'zoom',''};
  hpos = [240:-20:60];
  % get initial parameter values from st.vols{volhandle}.premul
  if volhandle == 0 
    volhandle = handles;
    prms(7:9) = 1;
  else
    prms = spm_imatrix(st.vols{volhandle}.premul);
    prms(10) = 3; % default #contour lines
    labels{end+1} = '#contour lines';
    st.vols{volhandle(1)}.reorient.b(1) = uicontrol(...
	Finter, 'Style','PushButton', 'Position',[75 30 165 025], ...
	'String','Apply to image(s)', ...
	'Callback',['spm_orthviews(''reorient'',''apply'',',...
		    num2str(volhandle), ');']);
  end;
  for k = handles
    obj = findobj(Fgraph, 'Tag',  ['REORIENT_M_', num2str(k)]);
    if any(k == volhandle)
      objh = findobj(obj, 'Tag', ['REORIENT_0_', num2str(k)]);
      objs = findobj(obj, 'Tag', ['REORIENT_1_', num2str(k)]);
      set(objh,'Visible','off');
      set(objs, 'Callback', ...
	       ['spm_orthviews(''reorient'',''context_quit'', [', ...
		num2str(volhandle), ']);'],'Visible','on');
      st.vols{k}.reorient.oldpremul = st.vols{k}.premul;
    else
      set(obj, 'Visible', 'off');
    end;
  end;
  for k = 1:numel(labels)
    st.vols{volhandle(1)}.reorient.l(k)=uicontrol(...
	Finter, 'Style','Text', ...
	'Position',[75 hpos(k) 100 016], 'String',labels{k});
    st.vols{volhandle(1)}.reorient.e(k) = uicontrol(...
	Finter, 'Style','edit', ...
	'Callback',['spm_orthviews(''reorient'',''reorient'',[',...
		    num2str(volhandle),'])'], ... 
	'Position',[175 hpos(k) 065 020], 'String',num2str(prms(k)), ...
	'ToolTipString',tooltips{k});
  end;
  spm_orthviews('redraw');
  
 case 'context_quit'
  Finter = spm_figure('FindWin', 'Interactive');
  Fgraph = spm_figure('FindWin', 'Graphics');
  try
    delete(st.vols{volhandle(1)}.reorient.e);
    delete(st.vols{volhandle(1)}.reorient.l);
    delete(st.vols{volhandle(1)}.reorient.b);
  catch
  end;
  if isfield(st.vols{volhandle(1)}.reorient,'lh')
    if ~isempty(st.vols{volhandle(1)}.reorient.lh)
      delete(cat(1,st.vols{volhandle(1)}.reorient.lh{:}));
    end;
  end;
  
  for k = spm_orthviews('valid_handles')
    try
      st.vols{k}.premul = st.vols{k}.reorient.oldpremul;
      st.vols{k} = rmfield(st.vols{k},'reorient');
    catch
    end;
    obj = findobj(Fgraph, 'Tag',  ['REORIENT_M_', num2str(k)]);
    if any(k == volhandle)
      objh = findobj(obj, 'Tag', ['REORIENT_1_', num2str(k)]);
      set(objh, 'Visible', 'off', 'Callback','');
      objs = findobj(obj, 'Tag', ['REORIENT_0_', num2str(k)]);
      set(objs, 'Visible', 'on');
    else
      set(obj, 'Visible', 'on');
    end;
  end;
  spm_orthviews('redraw');
  
  %-------------------------------------------------------------------------
  % Interaction callbacks
  
 case 'apply'
  [p n e v] = fileparts(st.vols{volhandle}.fname);
  P = cellstr(spm_select(Inf, 'image', {'Image(s) to reorient'}, '', p));
  if ~isempty(P)
    spm_progress_bar('Init', numel(P), 'Reorient', 'Images completed');
    for k = 1:numel(P)
      M = spm_get_space(P{k});
      spm_get_space(P{k},st.vols{volhandle}.premul*M);
      spm_progress_bar('Set',k);
    end;
    spm_progress_bar('Clear');
    st.vols{volhandle}.reorient.oldpremul = eye(4);
    qu=questdlg({'Image positions are changed!', ...
		 'To make sure images are displayed correctly, it is recommended to quit and restart spm_orthviews now.', ... 
		 'Do you want to quit?'},'Reorient done','Yes','No','Yes');
    if strcmp(lower(qu), 'yes')
      spm_orthviews('reset');
      return;
    end;
  end;
  spm_orthviews('reorient','context_quit', volhandle);
  
 case 'reorient'
  prms=zeros(1,12);
  for k=1:9
    prms(k) = str2num(get(st.vols{volhandle(1)}.reorient.e(k),'string'));
  end;
  for k = volhandle
    st.vols{k}.premul = spm_matrix(prms)* ...
	st.vols{k}.reorient.oldpremul;
  end;
  spm_orthviews('redraw');
  
 case 'redraw'
  if numel(st.vols{volhandle}.reorient.e)==10
    if isfield(st.vols{volhandle}.reorient,'lh')
      if ~isempty(st.vols{volhandle}.reorient.lh)
	delete(cat(1,st.vols{volhandle}.reorient.lh{:}));
      end;
    end;
    st.vols{volhandle}.reorient.lh = {};
    ncl = str2num(get(st.vols{volhandle}.reorient.e(10),'string'));
    if ncl > 0
      todraw=spm_orthviews('valid_handles');
      for d = 1:3
	CData = sqrt(sum(get(st.vols{volhandle}.ax{d}.d,'CData').^2, ...
			 3));
	for h = todraw
	  if h ~= volhandle
	    axes(st.vols{h}.ax{d}.ax);
	    hold on;
	    [C st.vols{volhandle}.reorient.lh{end+1}]=contour(CData,ncl,'r-');
	  end;
	end;
      end;
      set(cat(1,st.vols{volhandle}.reorient.lh{:}),'HitTest','off');
    end;
  end;
 otherwise    
  fprintf('spm_orthviews(''reorient'', ...): Unknown action %s', cmd);
end;

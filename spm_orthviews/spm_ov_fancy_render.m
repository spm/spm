function ret = spm_ov_fancy_render(varargin)
% Wrapper for tbxrend_fancy_render
% This functions provides an interface to tbxrend_fancy_render from
% an spm_orthviews display. The main features include:
% - initialisation based on the current spm_orthviews display
%   (volume data and blobs)
% - common crosshairs for better orientation
% - 3 high level cut/slice operations that cut out parts of the
%   surface from the current crosshair position towards the current
%   camera position:
%   - Cut+Slice box: Cut out a box where Cameraposition and current
%     crosshair position are opposite corners
%   - Cut+Slice y-interval box: crosshair y position is the middle of an
%     interval along the y axis, x and z positions on the edge
%   - Cut+Slice sphere: crosshairs on the surface of a sphere with
%     a pre-specified radius.
%
% This routine is a plugin to spm_orthviews for SPM5. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the matlab prompt.
%_______________________________________________________________________
%
% @(#) $Id: spm_ov_fancy_render.m,v 1.11 2006/02/23 10:34:28 glauche Exp $

rev = '$Revision: 1.11 $';

global st;
if isempty(st)
  error('fancy_render: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
  error('fancy_render: Wrong number of arguments. Usage: spm_orthviews(''fancy_render'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};

switch cmd
  
  %-------------------------------------------------------------------------
  % Context menu and callbacks
case 'context_menu'
    if exist(fullfile(spm('dir'),'toolbox','Render3D'),'dir')
        addpath(fullfile(spm('dir'),'toolbox','Render3D'));
    else
        return;
    end;
  if ~any(exist('tbxrend_fancy_render')==[2:6])
    warning([mfilename ':init'],['function tbxrend_fancy_render not' ...
		    ' found!']);
    return;
  end;
  VMfancy = tbxrend_fancy_render('menu',varargin{3});
  item0 = uimenu(VMfancy, 'Label','Cut+Slice box', 'Callback', ...
		 ['feval(''' mfilename ''',''cut_box'',' num2str(volhandle) ');'], ...
		 'Tag','VMFANCY1', 'Visible','off');
  item0 = uimenu(VMfancy, 'Label','Cut+Slice y interval', 'Callback', ...
		 ['feval(''' mfilename ''',''cut_ybox'',' num2str(volhandle) ');'], ...
		 'Tag','VMFANCY1', 'Visible','off');
  item0 = uimenu(VMfancy, 'Label','Cut+Slice sphere', 'Callback', ...
		 ['feval(''' mfilename ''',''cut_sphere'',' num2str(volhandle) ');'], ...
		 'Tag', 'VMFANCY1', 'Visible','off');
  VMfancyxh = uimenu(VMfancy, 'Label', 'Xhairs', 'Tag','VMFANCY1', ...
		     'Visible','off');
  VMfxhon = uimenu(VMfancyxh, 'Label','On', 'Checked','on',...
		   'Callback',...
		   ['feval(''' mfilename ''',''menuxhairson'',' num2str(volhandle) ');' ],...
		   'Tag','VMFANCY1_XHAIRSON');
  VMfxhoff = uimenu(VMfancyxh, 'Label','Off',...
		    'Callback',...
		    ['feval(''' mfilename ''',''menuxhairsoff'',' num2str(volhandle) ');' ],...
			    'Tag','VMFANCY1_XHAIRSOFF');    
  VMfancy0 = findobj(VMfancy, 'Tag','VMFANCY0');
  for k = 1:numel(VMfancy0)
    cb = get(VMfancy0(k),'Callback');
    cb = strrep(cb,'tbxrend_fancy_render',  mfilename);
    cb = strrep(cb,')',[',' num2str(volhandle) ')']);
    set(VMfancy0(k),'Callback',cb);    
  end;
  item1 = uimenu(VMfancy, 'Label', 'Help', 'Callback', ...
	  ['feval(''spm_help'',''' mfilename ''');']);
  
 case 'menuinit'
  fr = tbxrend_fancy_render('defaults');
  if isfield(st.vols{volhandle},'blobs')
    fr.blobs = st.vols{volhandle}.blobs;
    for k = 1:numel(fr.blobs)
            fr.blobs{k}.interp = 1;
            fr.blobs{k}.scaling = 'abs';
    end;
  end;
  fr.VV = spm_vol(st.vols{volhandle}.fname);
  fr = tbxrend_fancy_render('menuinit', fr);
  spm_ovhelper_3Dreg('register', 'tbxrend_fancy_render',fr.VV);
 
 case 'menuload'
  fr = tbxrend_fancy_render('menuload');
  spm_ovhelper_3Dreg('register', 'tbxrend_fancy_render',fr.VV);
 
 case 'menuxhairson'
  spm('pointer','watch');
  if nargin == 1
    mparent = get(gcbo,'Parent');
  else
    mparent = varargin{2};
  end;
  set(findobj(mparent, 'Tag','VMFANCY1_XHAIRSON'),'Checked','on');
  set(findobj(mparent, 'Tag','VMFANCY1_XHAIRSOFF'),'Checked','off');
  spm_ovhelper_3Dreg(varargin{1}(5:end), 'tbxrend_fancy_render');
  spm('pointer','arrow');
  return;
 
 case 'menuxhairsoff'
  spm('pointer','watch');
  if nargin == 1
    mparent = get(gcbo,'Parent');
  else
    mparent = varargin{2};
  end;
  set(findobj(mparent,'Tag','VMFANCY1_XHAIRSON'),'Checked','off');
  set(findobj(mparent,'Tag','VMFANCY1_XHAIRSOFF'),'Checked','on');  
  spm_ovhelper_3Dreg(varargin{1}(5:end), 'tbxrend_fancy_render');
  spm('pointer','arrow');
  return;
 
 case 'cut_box'
  pos = spm_orthviews('pos')';
  ax = findobj(0, 'Tag', 'tbxrend_fancy_render');
  fr = get(ax, 'Userdata');
  cp = get(ax,'Cameraposition');
  ct = get(ax, 'Cameratarget');
  cv = (cp-ct)./sqrt(sum((cp-ct).^2));
  cm = '>>>';
  cm(sign(cv)~=1) = '<';
  d  = 'xyz';
  ind = numel(fr.slicescuts)+1;
  fr.slicescuts{ind}.type = 'cut';
  for k = 1:3
    bounds{k} = sprintf('(%s %s= %d)', d(k), cm(k), pos(k));
  end;
  fr.slicescuts{ind}.val = sprintf('%s & %s & %s',...
				   bounds{:});
  for k = 1:3
    fr.slicescuts{ind+k}.type = 'slice';
    fr.slicescuts{ind+k}.val = sprintf('%s - %d', d(k), pos(k));
    fr.slicescuts{ind+k}.bounds = sprintf('%s & %s', bounds{[1 2 3]~=k});
  end;
  set(ax,'Userdata',fr);
  tbxrend_fancy_render('redraw');
 
 case 'cut_ybox'
  pos = spm_orthviews('pos')';
  yi = spm_input('Y interval (mm)','+1','e',[pos(2)-10 pos(2)+10],2);
  ax = findobj(0, 'Tag', 'tbxrend_fancy_render');
  fr = get(ax, 'Userdata');
  cp = get(ax,'Cameraposition');
  ct = get(ax, 'Cameratarget');
  cv = (cp-ct)./sqrt(sum((cp-ct).^2));
  cm = '>>>';
  cm(sign(cv)~=1) = '<';
  d  = 'xyz';
  ind = numel(fr.slicescuts)+1;
  fr.slicescuts{ind}.type = 'cut';
  for k = [1 3]
    bounds{k} = sprintf('(%s %s= %d)', d(k), cm(k), pos(k));
  end;
  bounds{2} = sprintf('(y >= %d) & (y <= %d)', yi);
  fr.slicescuts{ind}.val = sprintf('%s & %s & %s',...
				   bounds{:});
  for k = 1:4
    fr.slicescuts{ind+k}.type = 'slice';
  end;
  fr.slicescuts{ind+1}.val = sprintf('%s - %d', d(1), pos(1));
  fr.slicescuts{ind+1}.bounds = sprintf('%s & %s', bounds{[2 3]});
  fr.slicescuts{ind+2}.val = sprintf('%s - %d', d(2), yi(1));
  fr.slicescuts{ind+2}.bounds = sprintf('%s & %s', bounds{[1 3]});
  fr.slicescuts{ind+3}.val = sprintf('%s - %d', d(2), yi(2));
  fr.slicescuts{ind+3}.bounds = sprintf('%s & %s', bounds{[1 3]});
  fr.slicescuts{ind+4}.val = sprintf('%s - %d', d(3), pos(3));
  fr.slicescuts{ind+4}.bounds = sprintf('%s & %s', bounds{[1 2]});
  set(ax,'Userdata',fr);
  tbxrend_fancy_render('redraw');
 
 case 'cut_sphere'
  pos = spm_orthviews('pos')';
  r = spm_input('Radius of sphere (mm)','+1','e');
  ax = findobj(0, 'Tag', 'tbxrend_fancy_render');
  fr = get(ax, 'Userdata');
  cp = get(ax,'Cameraposition');
  ct = get(ax, 'Cameratarget');
  cv = (cp-ct)./sqrt(sum((cp-ct).^2));
  ind = numel(fr.slicescuts)+1;
  fr.slicescuts{ind}.type = 'cut';
  cent=r*cv+pos;
  sp = sprintf('(x-(%d)).^2+(y-(%d)).^2+(z-(%d)).^2 - (%d).^2', cent,r);
  fr.slicescuts{ind}.val = sprintf('%s<0',sp);
  fr.slicescuts{ind+1}.type = 'slice';
  fr.slicescuts{ind+1}.val = sp;
  set(ax,'Userdata',fr);
  tbxrend_fancy_render('redraw');
 
 otherwise    
  fprintf('spm_orthviews(''fancy_render'', ...): Unknown action %s', cmd);
end;

spm('pointer','arrow');

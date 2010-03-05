function spm_ovhelper_3Dreg(cmd, varargin)
% Helper function to register spm_orthviews plugins via spm_XYZreg
% FORMAT spm_ovhelper_3Dreg('register', h, V)
% Register a (3D) graphics with the main spm_orthviews display. This will
% draw 3D crosshairs at the current spm_orthviews position and update
% them whenever the spm_orthviews cursor moves.
% h - a graphics handle or a tag of graphics handle to register
% V - a volume handle (or equivalent) containing dimensions and
%     voxel-to-world mapping information
% FORMAT spm_ovhelper_3Dreg('unregister', h)
% h - a graphics handle or a tag of graphics handle to unregister
% FORMAT spm_ovhelper_3Dreg('setcoords', xyz, h)
% Update position of crosshairs in 3D display
% xyz - new crosshair coordinates (in mm)
% h - a graphics handle or a tag of graphics handle to update
% FORMAT spm_ovhelper_3Dreg('xhairson', h)
% FORMAT spm_ovhelper_3Dreg('xhairsoff', h)
% Toggle display of crosshairs in 3D display.
% h - a graphics handle or a tag of graphics handle 
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_ovhelper_3Dreg.m 3756 2010-03-05 18:43:37Z guillaume $

if ishandle(varargin{1})
  h = varargin{1};
elseif ischar(varargin{1})
  h = findobj(0, 'Tag',varargin{1});
  if ~ishandle(h)
    warning([mfilename ':InvalidHandle'], ...
        'No valid graphics handle found');
    return;
  else
    h = get(h(ishandle(h)),'parent');
  end;
end;

switch lower(cmd)
 case 'register'
  register(h,varargin{2:end});
  return;
 case 'setcoords'
  setcoords(varargin{1:end});
  return;
 case 'unregister',
  unregister(h,varargin{2:end});
  return;
 case 'xhairson'
  xhairs(h,'on',varargin{2:end});
  return;
 case 'xhairsoff'
  xhairs(h,'off',varargin{2:end});
  return;
end;

function register(h,V,varargin)
try
  global st;
  if isstruct(st)
    xyz=spm_orthviews('pos');
    if isfield(st,'registry')
      hreg = st.registry.hReg;
    else
      [hreg xyz]=spm_XYZreg('InitReg', h, V.mat, ...
                            V.dim(1:3)',xyz); 
      spm_orthviews('register',hreg);
    end;
    spm_XYZreg('Add2Reg',hreg,h,mfilename);
    feval(mfilename,'setcoords',xyz,h);
    set(h, 'DeleteFcn', ...
       sprintf('%s(''unregister'',%f);', mfilename, h));
  end;
catch
  warning([mfilename ':XYZreg'],...
          'Unable to register to spm_orthviews display');
  disp(lasterror);
end;    
return;

function setcoords(xyz,h,varargin)
spm('pointer','watch');
Xh = findobj(h,'Tag', 'Xhairs');
if ishandle(Xh)
  vis = get(Xh(1),'Visible');
  delete(Xh);
else
  vis = 'on';
end;
axes(findobj(h,'Type','axes'));
lim = axis;
Xh = line([xyz(1), xyz(1), lim(1);...
       xyz(1), xyz(1), lim(2)],...
      [lim(3), xyz(2), xyz(2);...
       lim(4), xyz(2), xyz(2)],...
      [xyz(3), lim(5), xyz(3);...
       xyz(3), lim(6), xyz(3)],...
      'Color','b', 'Tag','Xhairs', 'Visible',vis,...
      'Linewidth',2, 'HitTest','off');
spm('pointer','arrow');
return;

function xhairs(h,val,varargin)
Xh = findobj(h, 'Tag', 'Xhairs');
if ~isempty(Xh)
  set(Xh,'Visible',val);
end;
  
function unregister(h,varargin)
try
  global st;
  if isfield(st,'registry')
    hreg = st.registry.hReg;
  else
    hreg = findobj(0,'Tag','hReg');
  end;
  if h == hreg
    spm_XYZreg('UnInitReg',hreg);
    st = rmfield(st, 'registry');
  else
    spm_XYZreg('Del2Reg',hreg,h);
  end;
catch
  warning([mfilename ':XYZreg'],...
      'Unable to unregister');
  disp(lasterror);
end;
return;

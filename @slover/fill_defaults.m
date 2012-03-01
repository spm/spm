function obj = fill_defaults(obj)
% check and fill fields in object
% FORMAT obj = fill_defaults(obj)
% 
% Input
% obj    - object to fill
% 
% Output
% obj    - object filled
%
% $Id: fill_defaults.m,v 1.2 2005/05/06 22:59:56 matthewbrett Exp $

% Some default structures
def_labs = struct('colour',[1 1 1],'size',0.075,'format', '%+3.0f');
def_fig  = struct('position', [0 0 1 0.92], 'units', 'normalized', ...
             'valign', 'top');
def_area = struct('position', [0 0 1 1], ...
          'units', '', ...
          'halign', 'center',...
          'valign', 'middle');

% Figure.  We allow the figure to be dead, if we are going to resurrect
% it later.
dead_f = 0;
if ~isempty(obj.figure)
  % Is it dead?
  if ~ishandle(obj.figure)
    % Do we want to revive it?
    if ~obj.resurrectf
      error('Figure handle is not a valid figure')
    end
    dead_f = 1;
    obj.refreshf = 1;
  elseif ~strcmp(get(obj.figure,'Type'),'figure')
    error('Figure handle is not a figure')
  end
else
  % no figure handle. Try spm figure, then gcf
  obj.figure = spm_figure('FindWin', 'Graphics'); 
  if isempty(obj.figure)
    obj.figure = gcf;
  end
end
% set defaults for SPM figure 
if ~dead_f
  if strcmp(get(obj.figure, 'Tag'),'Graphics')
    % position figure nicely for SPM
    obj.area = mars_struct('fillafromb', obj.area, def_fig);
  end
end

% orientation; string or 4x4 matrix
if ischar(obj.transform)
  orientn = find(strcmpi(obj.transform, {'axial', ... 
            'coronal', ...
            'sagittal'}));
  if isempty(orientn)
    error('Unexpected orientation %s', obj.transform);
  end
  ts = [0 0 0 0 0 0 1 1 1;...
      0 0 0 pi/2 0 0 1 -1 1;...
      0 0 0 pi/2 0 -pi/2 -1 1 1];
  obj.transform = spm_matrix(ts(orientn,:));
end

% default slice size, slice matrix depends on orientation
if (isempty(obj.slicedef) || isempty(obj.slices)) ...
      && ~isempty(obj.img)
  % take image sizes from first image
  V = obj.img(1).vol;
  D = V.dim(1:3);
  T = obj.transform * V.mat;
  vcorners = [1 1 1; D(1) 1 1; 1 D(2) 1; D(1:2) 1; ...
          1 1 D(3); D(1) 1 D(3); 1 D(2:3) ; D(1:3)]';
  corners = T * [vcorners; ones(1,8)];
  SC = sort(corners, 2);
  vxsz = sqrt(sum(T(1:3,1:3).^2));
  
  if isempty(obj.slicedef)
    obj.slicedef = [SC(1,1) vxsz(1) SC(1,8);SC(2,1) vxsz(2) SC(2,8)];
  end
  if isempty(obj.slices)
    obj.slices = SC(3,1):vxsz(3):SC(3,8);
  end
end

% labels
if ischar(obj.labels)
  if ~strcmpi(obj.labels, 'none')
    error('If labels is string, should be ''none''');
  end
else
  obj.labels = mars_struct('fillafromb', obj.labels, def_labs);
end

% figure area stuff
obj.area = mars_struct('fillafromb', obj.area, def_area);
if isempty(obj.area.units)
  if (all(obj.area.position>=0 & obj.area.position<=1))
    obj.area.units = 'normalized';
  else
    obj.area.units = 'pixels';
  end
end

% fill various img arguments

% set colour intensities as we go
remcol = 1;
for i = 1:length(obj.img)
  if ~mars_struct('isthere', obj.img(i), 'type')
    % default is true colour, unless prop is Inf
    obj.img(i).type = 'truecolour';
    if mars_struct('isthere', obj.img(i), 'prop')
      if obj.img(i).prop == Inf
    obj.img(i).type = 'split';
    obj.img(i).prop = 1;
      end
    end
  end
  if ~mars_struct('isthere', obj.img(i), 'hold')
    if ~mars_struct('isthere', obj.img(i).vol, 'imgdata')
      % normal file vol struct
      obj.img(i).hold = 1;
    else
      % 3d matrix vol struct
      obj.img(i).hold = 0;
    end
  end
  if ~mars_struct('isthere', obj.img(i), 'background')
    obj.img(i).background = NaN;
  end
  if ~mars_struct('isthere', obj.img(i), 'prop')
    % default is true colour
    if strcmpi(obj.img(i).type, 'truecolour')
      obj.img(i).prop = remcol/(length(obj.img)-i+1);
      remcol = remcol - obj.img(i).prop;
    else
      obj.img(i).prop = 1;
    end
  end
  if ~mars_struct('isthere', obj.img(i), 'range')
    [mx mn] = pr_volmaxmin(obj.img(i).vol);
    obj.img(i).range = [mn mx];
  end
  if ~mars_struct('isthere', obj.img(i), 'cmap')
    if strcmpi(obj.img(i).type, 'split')
      if obj.range(1)<obj.range(2)
    obj.img(i).cmap = pr_getcmap('hot');
      else
    obj.img(i).cmap = pr_getcmap('winter');
      end
    else                  % true colour
      obj.img(i).cmap = gray;
    end
  else % check cmap is OK
    if ischar(obj.img(i).cmap)
      obj.img(i).cmap = pr_getcmap(obj.img(i).cmap);
    end
  end  
  if ~mars_struct('isthere', obj.img(i), 'outofrange')
    % this can be complex, and depends on split/true colour
    if strcmpi(obj.img(i).type, 'split')
      if xor(obj.img(i).range(1) < obj.img(i).range(2), ...
         obj.img(i).range(2) < 0)
    obj.img(i).outofrange = {0, size(obj.img(i).cmap,1)};
      else
    obj.img(i).outofrange={1, 0};
      end
    else            % true colour
      obj.img(i).outofrange = {1,size(obj.img(i).cmap,1)};
    end
  end
  for j=1:2
    if isempty(obj.img(i).outofrange{j})
      obj.img(i).outofrange(j) = {0};
    end
  end
  if ~mars_struct('isthere', obj.img(i), 'nancol')
    obj.img(i).nancol = 0;
  end
end
return

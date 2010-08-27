function [ans1, ans2, ans3] = axis(varargin)
%AXIS  Control axis scaling and appearance.
%   AXIS([XMIN XMAX YMIN YMAX]) sets scaling for the x- and y-axes
%      on the current plot.
%   AXIS([XMIN XMAX YMIN YMAX ZMIN ZMAX]) sets the scaling for the
%      x-, y- and z-axes on the current 3-D plot.
%   AXIS([XMIN XMAX YMIN YMAX ZMIN ZMAX CMIN CMAX]) sets the
%      scaling for the x-, y-, z-axes and color scaling limits on
%      the current axis (see CAXIS). 
%   V = AXIS returns a row vector containing the scaling for the
%      current plot.  If the current view is 2-D, V has four
%      components; if it is 3-D, V has six components.
%
%   AXIS AUTO  returns the axis scaling to its default, automatic
%      mode where, for each dimension, 'nice' limits are chosen based
%      on the extents of all line, surface, patch, and image children.
%   AXIS MANUAL  freezes the scaling at the current limits, so that if
%      HOLD is turned on, subsequent plots will use the same limits.
%   AXIS TIGHT  sets the axis limits to the range of the data.
%   AXIS FILL  sets the axis limits and PlotBoxAspectRatio so that
%      the axis fills the position rectangle.  This option only has
%      an effect if PlotBoxAspectRatioMode or DataAspectRatioMode are
%      manual.
%
%   AXIS IJ  puts MATLAB into its "matrix" axes mode.  The coordinate
%      system origin is at the upper left corner.  The i axis is
%      vertical and is numbered from top to bottom.  The j axis is
%      horizontal and is numbered from left to right.
%   AXIS XY  puts MATLAB into its default "Cartesian" axes mode.  The
%      coordinate system origin is at the lower left corner.  The x
%      axis is horizontal and is numbered from left to right.  The y
%      axis is vertical and is numbered from bottom to top.
%
%   AXIS EQUAL  sets the aspect ratio so that equal tick mark
%      increments on the x-,y- and z-axis are equal in size. This
%      makes SPHERE(25) look like a sphere, instead of an ellipsoid.
%   AXIS IMAGE  is the same as AXIS EQUAL except that the plot
%      box fits tightly around the data.
%   AXIS SQUARE  makes the current axis box square in size.
%   AXIS NORMAL  restores the current axis box to full size and
%       removes any restrictions on the scaling of the units.
%       This undoes the effects of AXIS SQUARE and AXIS EQUAL.
%   AXIS VIS3D  freezes aspect ratio properties to enable rotation of
%       3-D objects and overrides stretch-to-fill.
%
%   AXIS OFF  turns off all axis labeling, tick marks and background.
%   AXIS ON  turns axis labeling, tick marks and background back on.
%
%   AXIS(H,...) changes the axes handles listed in vector H.
%
%   See also AXES, GRID, SUBPLOT, XLIM, YLIM, ZLIM.

%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 4052 $  $Date: 2008/05/05 21:38:11 $

%get the list of axes to operate upon
if ~isempty(varargin) && allAxes(varargin{1})
    
    ax = varargin{1};
    varargin=varargin(2:end);
else
    ax = gca;
end

ans1set = false;
pbarlimit = 0.1;

%---Check for bypass option (only supported for single axes)
if length(ax)==1 && isappdata(ax,'MWBYPASS_axis')
    if isempty(varargin)
        ans1 = mwbypass(ax,'MWBYPASS_axis');
        ans1set = true;
    else
        mwbypass(ax,'MWBYPASS_axis',varargin{:});
    end
elseif isempty(varargin)
    if length(ax)==1
        ans1=LocGetLimits(ax);
        ans1set = true;
    else
        ans1=cell(length(ax),1);
        ans1set = true;
        for i=1:length(ax)
            ans1{i}=LocGetLimits(ax(i));    
        end
    end
else
    for j=1:length(ax)
        for i = 1:length(varargin)
            cur_arg = varargin{i};
            
            % Set limits manually with 4/6/8 element vector
            if isnumeric(cur_arg)
                LocSetLimits(ax(j),cur_arg);
                
                % handle AUTO, AUTO[XYZ]:
            elseif strcmp(cur_arg(1:min(4,length(cur_arg))),'auto')
                LocSetAuto(ax(j),cur_arg);
                
                % handle TIGHT
            elseif(strcmp(cur_arg,'tight'))
                LocSetTight(ax(j));
                
                % handle FILL:
            elseif(strcmp(cur_arg, 'fill'))
                LocSetFill(ax(j),pbarlimit);
                
                % handle MANUAL:
            elseif(strcmp(cur_arg, 'manual'))
                set(ax(j),...
                    'XLimMode','manual',...
                    'YLimMode','manual',...
                    'ZLimMode','manual');
                
                % handle IJ:
            elseif(strcmp(cur_arg, 'ij'))
                set(ax(j),...
                    'XDir','normal',...
                    'YDir','reverse');
                
                % handle XY:
            elseif(strcmp(cur_arg, 'xy'))
                set(ax(j),...
                    'XDir','normal',...
                    'YDir','normal');
                
                % handle SQUARE:
            elseif(strcmp(cur_arg, 'square')) 
                set(ax(j),...
                    'PlotBoxAspectRatio',[1 1 1],...
                    'DataAspectRatioMode','auto')
                
                % handle EQUAL:
            elseif(strcmp(cur_arg, 'equal'))
                LocSetEqual(ax(j),pbarlimit);
                
                % handle IMAGE:
            elseif(strcmp(cur_arg,'image'))
                LocSetImage(ax(j),pbarlimit);
                
                % handle NORMAL:
            elseif(strcmp(cur_arg, 'normal'))
                set(ax(j),...
                    'PlotBoxAspectRatioMode','auto', ...
                    'DataAspectRatioMode'   ,'auto', ...
                    'CameraViewAngleMode'   ,'auto');
                
                % handle OFF:
            elseif(strcmp(cur_arg, 'off'))
                set(ax(j),'Visible','off');
                set(get(ax(j),'Title'),'Visible','on');
                
            % handle ON:
            elseif(strcmp(cur_arg, 'on'))
                set(ax(j),'Visible','on');
                
            % handle VIS3D:
            elseif(strcmp(cur_arg,'vis3d'))
                set(ax(j),...
                    'CameraViewAngle',   get(ax(j),'CameraViewAngle'),...
                    'DataAspectRatio',   get(ax(j),'DataAspectRatio'),...
                    'PlotBoxAspectRatio',get(ax(j),'PlotBoxAspectRatio'));
                
            % handle STATE:
            elseif(strcmp(cur_arg, 'state'))
                warning('MATLAB:graph2d:axis:ObsoleteState',(['AXIS(''STATE'') is obsolete and', ...
                        ' will be eliminated\n         in future versions.', ...
                        ' Use GET(GCA,...) instead.']));
                %note that this will keep overwriting arg1 etc if there is more
                %than one axes in the list
                [ans1,ans2,ans3]=LocGetState(ax(1));
                ans1set = true;
                
                %if nargout>1
                %    ans2=ans2q;
                %    if nargout>2
                %        ans3=ans3q;
                %    end
                %end
                
            % handle ERROR (NONE OF THE ABOVE STRINGS FOUND):
            else
                error('MATLAB:axis:UnknownOption', 'Unknown command option %s', cur_arg);
            end
        end
    end
end

if nargout > 0 && ~ans1set
    error(nargoutchk(0, 0, nargout));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ans1=LocGetLimits(axH)
%returns a 4 or 6 element vector of limits for a single axis

ans1 = [get(axH,'XLim') get(axH,'YLim')];
v = get(axH,'View');
if (v(1) ~= 0 && v(2) ~= 90)  % it's effectively 2D (v ~= [0 90])
    ans1 = [ans1 get(axH,'ZLim')];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocSetLimits(ax,lims)

if (length(lims) == 4) || (length(lims) == 6 || (length(lims) == 8))
    set(ax,...
        'XLim',lims(1:2),...
        'YLim',lims(3:4),...
        'XLimMode','manual',...
        'YLimMode','manual');
    
    if length(lims) == 6 || length(lims) == 8
        set(ax,...
            'ZLim',lims(5:6),...
            'ZLimMode','manual');
    end
    
    if length(lims) == 8
        set(ax,...
            'CLim',lims(7:8),...
            'CLimMode','manual');
    end
    
    if length(lims) == 4 && ~ishold(ax)
        set(ax,'view',[0 90]);
    elseif (length(lims) == 6 && ...
            isequal(get(ax,'View'),[0 90]) && ...
            ~ishold(ax))
        set(ax,'view',[-37.5,30]);
    end
else
    error('MATLAB:axis:WrongNumberElements', 'Vector must have 4, 6, or 8 elements.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocSetAuto(ax,cur_arg)
%called in response to axis auto[xyz]

do_all = (length(cur_arg) == length('auto'));
do_x = length(find(cur_arg == 'x'));
do_y = length(find(cur_arg == 'y'));
do_z = length(find(cur_arg == 'z'));
if(do_all || do_x)
    set(ax,'XLimMode','auto');
else
    set(ax,'XLimMode','manual');
end
if(do_all || do_y)
    set(ax,'YLimMode','auto');
else
    set(ax,'YLimMode','manual');
end
if(do_all || do_z)
    set(ax,'ZLimMode','auto');
else
    set(ax,'ZLimMode','manual');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocSetTight(ax)
%called in response to axis tight
%we can assume that ax is length=1

limits = objbounds(findall(ax));
if isempty(limits) % no objects in axes with data limits
  limits = [get(ax,'XLim') get(ax,'YLim') get(ax,'ZLim')];
end
hasdepth = limits(5) ~= limits(6);

% Protect against axis limit values being the same
ndx = find(diff(limits)==0 & [1 0 1 0 1]);

% handle log scales
logscales = [0 0 0 0 0 0];
logscales(1) = strcmp(get(ax,'xscale'),'log');
logscales(3) = strcmp(get(ax,'yscale'),'log');
if hasdepth
    logscales(5) = strcmp(get(ax,'zscale'),'log');
end

if ~isempty(ndx)
    if any(logscales(ndx))
        for i = 1:length(ndx)
            j = ndx(i);
            if logscales(i)
                % handle semilogx(1,1:10)
                % Scale upper 10 (log)
                % Scale lower limit by .1 (log)
                limits(j+1) = limits(j) * 10;
                limits(j) = limits(j) / 10;
            else
                % handle semilogy(1:10,1)
                % Bump upper and lower limits by 1
                limits(j+1) = limits(j)+1;
                limits(j) = limits(j)-1;
            end
        end
    else
        % Bump upper and lower limits by 1
        % handle semilogx(1:10,1)
        % handle semilogy(1,1:10)
        % handle plot(1,1:10), plot(1,1), plot(1:10,1)
        limits(ndx+1) = limits(ndx)+1;
        limits(ndx) = limits(ndx)-1;
    end
end

if all(isfinite(limits(1:4))),
    set(ax,...
        'xlim',limits(1:2),...
        'ylim',limits(3:4) + [-1 1]*diff(limits(3:4))/16)
end

if hasdepth
    set(ax,'zlim',limits(5:6))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocSetFill(ax,pbarlimit)
%called in response to axis fill

if strcmp(get(ax,'PlotBoxAspectRatioMode'),'manual') || ...
        strcmp(get(ax,'DataAspectRatioMode'),'manual')
    % Check for 3-D plot
    if all(rem(get(ax,'view'),90)~=0),
        a = axis;
        axis auto
        axis image
        pbar = get(ax,'PlotBoxAspectRatio');
        
        if pbar(1)~=pbarlimit, set(ax,'xlim',a(1:2)); end
        if pbar(2)~=pbarlimit, set(ax,'ylim',a(3:4)); end
        if pbar(3)~=pbarlimit, set(ax,'zlim',a(5:6)); end
        return
    end
    
    units = get(ax,'Units'); set(ax,'Units','Pixels')
    a = get(ax,'Position'); set(ax,'Units',units)
    % Change the unconstrained axis limit to 'auto'
    % based on the axis position.  Also set the pbar.
    set(ax,'PlotBoxAspectRatio',a([3 4 4]))
    if a(3) > a(4), 
        set(ax,'xlimmode','auto')
    else
        set(ax,'ylimmode','auto')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocSetEqual(ax,pbarlimit)
%called in response to axis equal

% Check for 3-D plot.  If so, use AXIS IMAGE.
if all(rem(get(ax,'view'),90)~=0),
    LocSetImage(ax,pbarlimit);
    return
end

units = get(ax,'Units'); set(ax,'Units','Pixels')
a = get(ax,'Position'); set(ax,'Units',units)
set(ax,'DataAspectRatio',[1 1 1]);
dx = diff(get(ax,'xlim')); dy = diff(get(ax,'ylim'));
dz = diff(get(ax,'zlim'));
set(ax,'PlotBoxAspectRatioMode','auto')
pbar = get(ax,'PlotBoxAspectRatio');
set(ax,'PlotBoxAspectRatio', ...
    [a(3) a(4) dz*min(a(3),a(4))/min(dx,dy)]);

% Change the unconstrained axis limit to auto based
% on the PBAR.
if pbar(1)/a(3) < pbar(2)/a(4), 
    set(ax,'xlimmode','auto')
else
    set(ax,'ylimmode','auto')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocSetImage(ax,pbarlimit)

set(ax,...
    'DataAspectRatio',[1 1 1], ...
    'PlotBoxAspectRatioMode','auto')

% Limit plotbox aspect ratio to 1 to 25 ratio.
pbar = get(ax,'PlotBoxAspectRatio');
pbar = max(pbarlimit,pbar / max(pbar));
if any(pbar(1:2) == pbarlimit),
    set(ax,'PlotBoxAspectRatio',pbar)
end

LocSetTight(ax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ans1,ans2,ans3]= LocGetState(ax)

str = '';
if(strcmp(get(ax,'XLimMode'), 'auto'))
    str = 'x';
end

if(strcmp(get(ax,'YLimMode'), 'auto'))
    str = [str, 'y'];
end

if(strcmp(get(ax,'ZLimMode'), 'auto'))
    str = [str, 'z'];
end

if length(str) == 3
    ans1 = 'auto';
else
    ans1 = 'manual';
end

if strcmp(get(ax,'Visible'),'on')
    ans2 = 'on';
else
    ans2 = 'off';
end

if(strcmp(get(ax,'XDir'),'normal') && ...
        strcmp(get(ax,'YDir'),'reverse'))
    ans3 = 'ij';
else
    ans3 = 'xy';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = allAxes(h)

result = all(ishghandle(h)) && ...
         length(findobj(h,'type','axes','-depth',0)) == length(h);

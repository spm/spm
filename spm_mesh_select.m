function P = spm_mesh_select(M,N)
% Select vertices interactively on a triangle mesh
% FORMAT P = spm_mesh_select(M,N)
% M        - a mesh filename or GIfTI object or patch structure
% N        - number of points to be interactively selected [default: 3]
%            or cell array of char vectors containing label of points
%
% P        - array of selected vertices coordinates [3xN]
%__________________________________________________________________________

% Tim Tierney, Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


%-Get input mesh
M = gifti(M);
if ~isfield(M,'cdata')
    M.cdata = repmat([0.7 0.7 0.7],size(M.vertices,1),1);
end

%-Get number of points to be selected
if nargin < 2
    N = 3;
end
if isnumeric(N)
    labels = arrayfun(@(i) sprintf('P%d',i),1:N,'UniformOutput',false);
else
    labels = N;
    N = numel(labels);
end

%-Display mesh in a new figure
F = figure;
H = spm_mesh_render(M,'Parent',F);
cameratoolbar(H.figure,'Show');
cameratoolbar(H.figure,'SetMode','orbit');
cameratoolbar(H.figure,'SetCoordSys','none');
set(H.light,'Visible','off');
camlight(H.axis,0,-90);
camlight(H.axis,-90,90);
hold(H.axis,'on');

%-Display buttons
for i=1:N
    H.button(i) = uicontrol(H.figure,...
        'Style','togglebutton',...
        'String',labels{i},...
        'Units','pixels',...
        'Position',[20,20,60,20].*[1 1+2*(i-1) 1 1],...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'Callback',{@cb_fig,i});

    H.box(i) = uicontrol(H.figure,...
        'Style','text',...
        'String','',...
        'Units','pixels',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'Position',[100,20,100,20].*[1 1+2*(i-1) 1 1]);

    H.dot(i) = plot3(H.axis,0,0,0,'r.','MarkerSize',20,'Visible','off');

end
setappdata(H.axis,'handles',H);

%-Global storage of selected vertices
setappdata(0,'vertex',NaN(3,N));
 
%-Wait until all vertices have been selected
waitfor(H.figure);

P = getappdata(0,'vertex');
rmappdata(0,'vertex');


%==========================================================================
function cb_fig(src, evt, idx)
if get(gcbo,'Value') == 1
    set(gcbf, 'WindowButtonDownFcn', {@cb_select_vertex, idx}); 
else
    set(gcbf, 'WindowButtonDownFcn', []);
end


%==========================================================================
function cb_select_vertex(src, evt, idx)

ax = ancestor(findobj(src,'Tag','SPMMeshRender'),'axes');
H  = getappdata(ax,'handles');

%-Get mesh vertices
vertices  = get(H.patch, 'Vertices')';

%-Get mouse selection type (e.g. to filter for double-click only)
selType   = get(src,'SelectionType');
if ~ismember(selType,{'normal'}) % switch to 'open' for double-click
    return;
end

%-Get camera properties
% See https://uk.mathworks.com/help/matlab/creating_plots/defining-scenes-with-camera-graphics.html
currPoint = get(H.axis, 'CurrentPoint')';
camPos    = get(H.axis, 'CameraPosition');
camTarget = get(H.axis, 'CameraTarget');
camUp     = get(H.axis, 'CameraUpVector');

%-Define transformation to the viewing frame
camDir = camPos - camTarget;
zAxis  = camDir / norm(camDir);
upAxis = camUp / norm(camUp);
xAxis  = cross(upAxis, zAxis);
yAxis  = cross(zAxis, xAxis);
T      = [xAxis; yAxis; zAxis];

%-Transform data into the viewing frame
VFvertices  = T * vertices;
VFcurrPoint = T * currPoint;

%-Find the nearest vertex on the mesh (using heuristic)
v = VFvertices;
nnearest = zeros(1,8);
for i=1:numel(nnearest)
    nnearest(i) = dsearchn(v(1:2,:)', VFcurrPoint(1:2));
    v(:,nnearest(i)) = Inf;
end
[~,i] = max(VFvertices(3,nnearest));
vertex = vertices(:, nnearest(i));

%-Store vertex coordinates
P = getappdata(0,'vertex');
P(:,idx) = vertex;
setappdata(0,'vertex',P);

%-Display selected vextex, print its coordinates and reset corresponding button
set(H.dot(idx),'XData',vertex(1,:),'YData',vertex(2,:),'ZData',vertex(3,:),'Visible','on');
set(H.box(idx),'String',sprintf('%d %d %d',round(vertex)));
set(H.button(idx),'Value',0);

%-Return to orbit mode (and implicitly remove the figure callback)
cameratoolbar('SetMode','orbit');
cameratoolbar('SetCoordSys','none');

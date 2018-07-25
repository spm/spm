function varargout = plot(varargin)
% plot method for GIfTI objects
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: plot.m 7381 2018-07-25 10:27:54Z guillaume $

% if ishandle(varargin{1})
%     h = figure(varargin{1});
% else
%     h = figure;
%     %axis equal;
%     %axis off;
%     %camlight;
%     %camlight(-80,-10);
%     %lighting gouraud;
% end
% cameratoolbar;


cdata = [];
ax = [];
if nargin == 1
    this = varargin{1};
    h = gcf;
else
    if ishandle(varargin{1})
        ax = varargin{1};
        h = figure(get(ax,'parent'));
        this = varargin{2};
    else
        this = varargin{1};
        h = gcf;
        cdata = subsref(varargin{2},struct('type','.','subs','cdata'));
    end
    if nargin > 2
        indc = varargin{3};
    else
        indc = 1;
    end
end

if isempty(ax), ax = axes('Parent',h); end
axis(ax,'equal');
axis(ax,'off');
hp = patch(struct(...
    'vertices',  subsref(this,struct('type','.','subs','vertices')),...
    'faces',     subsref(this,struct('type','.','subs','faces'))),...
    'FaceColor', 'b',...
    'EdgeColor', 'none',...
    'Parent',ax);

if ~isempty(cdata)
    set(hp,'FaceVertexCData',cdata(:,indc), 'FaceColor','interp')
end

camlight(ax);
camlight(ax,-80,-10);
lighting(ax,'gouraud');
cameratoolbar(h);

if nargout
    varargout{1} = hp;
end

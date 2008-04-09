function layoutplot(cfg, data);

% LAYOUTPLOT makes a figure with the 2-D layout of the channel positions
% for topoplotting and the individual channel axes (i.e. width and height
% of the subfigures) for multiplotting. A correct 2-D layout is a
% prerequisite  for plotting the topographical distribution of the
% potential or field distribution, or for plotting timecourses in a
% topographical arrangement.
%
% This function uses the same configuration options as prepare_layout and
% as the topoplotting and multiplotting functions. The difference is that
% this function only plots the layout without any data, which facilitates
% the validation of your 2-D layout.
%
% Use as
%   layoutplot(cfg, data)
%
% There are several ways in which a 2-D layout can be made: it can be read
% directly from a *.lay file, it can be created based on 3-D electrode or
% gradiometer positions in the configuration or in the data, or it can be
% created based on the specification of an electrode of gradiometer file.
%
% You can specify either one of the following configuration options
%   cfg.layout      filename containg the 2-D layout
%   cfg.rotate      number, rotation around the z-axis in degrees (default = [], which means automatic)
%   cfg.projection  string, 'stereographic', 'ortographic', 'polar', 'gnomic' or 'inverse' (default = 'orthographic')
%   cfg.elec        structure with electrode positions, or
%   cfg.elecfile    filename containing electrode positions
%   cfg.grad        structure with gradiometer definition, or
%   cfg.gradfile    filename containing gradiometer definition
%
% Alternatively the layout can be constructed from either
%   data.elec     structure with electrode positions
%   data.grad     structure with gradiometer definition
%
% See also prepare_layout, topoplotER, topoplotTFR, multiplotER, multiplotTFR

% Copyright (C) 2006-2007, Robert Oostenveld
%
% $Log: layoutplot.m,v $
% Revision 1.6  2007/03/21 14:47:06  chrhes
% updated some comments
%
% Revision 1.5  2007/03/20 10:44:37  roboos
% updated documentation, change figure axis, changed whitespaces
%
% Revision 1.4  2007/03/20 09:46:31  chrhes
% added some checks to determine whether the field cfg.layout already contains a
% valid layout structure that can be used, thereby avoiding a redundant call to
% the function PREPARE_LAYOUT
%
% Revision 1.3  2007/03/14 08:56:00  roboos
% changed some documentation
%
% Revision 1.2  2007/03/14 08:53:29  roboos
% changed from using private/createlayout to prepare_layout, changed the lay
% structure, added some help to clarify the usefullness of this function
%
% Revision 1.1  2006/06/01 11:51:45  roboos
% new implementation
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic check/initialization of input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<1) || (nargin>2), error('incorrect number of input arguments'); end;
if (nargin<2), data = []; end;

if ~isstruct(cfg), error('argument cfg must be a structure'); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract/generate layout information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lay = [];
% try to use the layout structure from input if specified
if isfield(cfg, 'layout')
  % brief check to determine if cfg.layout is a valid layout (lay) structre
  if isstruct(cfg.layout)
    if all(isfield(cfg.layout, {'pos';'width';'height';'label'}))
      lay = cfg.layout;
    end
  end
end
% otherwise create the layout structure
if isempty(lay), lay = prepare_layout(cfg, data); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot channel labels and local axes (boxes) in one figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X      = lay.pos(:,1);
Y      = lay.pos(:,2);
Width  = lay.width;
Height = lay.height;
Lbl    = lay.label;
figure
plot(X, Y, '.');
text(X, Y, Lbl);
line([X X+Width X+Width X X]',[Y Y Y+Height Y+Height Y]');
axis auto
axis equal

function [vol, sens] = headmodelplot(cfg, data);

% HEADMODELPLOT makes a 3D visualisation of the volume conductor model
% and optionally of the gradiometer positions and headshape. It can
% be used for example to check CTF multiple-sphere head models.
%
% Use as
%   headmodelplot(cfg)
%   headmodelplot(cfg, data)
%
% You should specify the volume conductor model with
%   cfg.hdmfile       = string, file containing the volume conduction model
% or alternatively
%   cfg.vol           = structure with volume conduction model
%
% If the sensor information is not contained in the data itself you should
% also specify the sensor information using
%   cfg.gradfile      = string, file containing the gradiometer definition
%   cfg.elecfile      = string, file containing the electrode definition
% or alternatively
%   cfg.grad          = structure with gradiometer definition
%   cfg.elec          = structure with electrode definition
%
% The positions of the sources can be specified as a regular 3-D
% grid that is aligned with the axes of the head coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
% Alternatively the position of a few sources at locations of interest can
% be specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = Nx3 matrix with position of each source
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%   cfg.grid.inside     = vector with indices of the sources inside the brain (optional)
%   cfg.grid.outside    = vector with indices of the sources outside the brain (optional)
% You can also use the PREPARE_LEADFIELD or SOURCEANALYSIS functions
% to create a grid with dipole positions.
%
% Other options are
%   cfg.channel          = cell-array, see CHANNELSELECTION
%   cfg.headshape        = name of CTF *.shape file
%   cfg.spheremesh       = number of vertices for spheres, either 42, 162 or 642
%   cfg.plotheadsurface  = 'yes' or 'no', is constructed from head model
%   cfg.plotheadshape    = 'yes' or 'no'
%   cfg.plotspheres      = 'yes' or 'no'
%   cfg.plotspherecenter = 'yes' or 'no'
%   cfg.plotgrid         = 'yes' or 'no'
%   cfg.plotinside       = 'yes' or 'no'
%   cfg.plotoutside      = 'yes' or 'no'
%   cfg.plotsensors      = 'yes' or 'no' plot electrodes or gradiometers
%   cfg.plotcoil         = 'yes' or 'no' plot all gradiometer coils
%   cfg.plotlines        = 'yes' or 'no' plot lines from sensor to head surface

% Undocumented local options:
% cfg.surface_facecolor
% cfg.surface_edgecolor
% cfg.surface_facealpha
%
% This function depends on PREPARE_VOL_SENS which has the following options:
% cfg.channel, documented
% cfg.elec, documented
% cfg.elecfile, documented
% cfg.grad, documented
% cfg.gradfile, documented
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, (default set in HEADMODELPLOT at cfg.vol =[]), documented

% Copyright (C) 2004-2007, Robert Oostenveld
%
% $Log: headmodelplot.m,v $
% Revision 1.19  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.18  2007/10/25 12:13:45  roboos
% make better determination of channel number for CTF gradiometer systems with broken channels
%
% Revision 1.17  2007/08/06 09:20:14  roboos
% added support for bti_hs
%
% Revision 1.16  2007/07/26 07:58:24  roboos
% also deal with cfg.headshape specified as surface, set of points or ctf_hs file (todo: add bti_hs)
%
% Revision 1.15  2007/07/26 07:07:08  roboos
% add empty inside and outside to empty sourcegrid
%
% Revision 1.14  2007/05/16 11:47:39  roboos
% added plotting of dipole grid, default for inside=yes and outside=no
% updated documentation
%
% Revision 1.13  2007/05/01 08:25:36  roboos
% fixed bug that was due to search-and-replace
%
% Revision 1.12  2006/10/04 08:19:42  roboos
% renamed megsystem into sensortype
%
% Revision 1.11  2006/07/24 07:59:16  roboos
% updated documentation
%
% Revision 1.10  2006/04/25 12:52:41  roboos
% added the facecolor and edgecolor options to the cfg as default
%
% Revision 1.9  2006/04/20 09:15:09  roboos
% updated documentation
% made color and opacity for the head surface externally available through cfg
%
% Revision 1.8  2006/04/10 16:33:46  ingnie
% updated documentation
%
% Revision 1.7  2005/12/14 10:44:23  roboos
% implemented support for EEG, spherical and BEM
% implemented support for single sphere in  case of MEG
% fixed some incompatibilities in grad and vol structure for MEG
% removed cfg.plotheadshape, just look whether isempty(cfg.headshape)
%
% Revision 1.6  2005/11/24 16:26:07  roboos
% updated multisphere selection to reflect the current change in multisphere vol structure (no labels any more)
%
% Revision 1.5  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.4  2004/08/19 07:22:49  roboos
% added cfg.plotspherecenter option
%
% Revision 1.3  2004/08/19 07:19:55  roboos
% added cfg.plotcoil option
%
% Revision 1.2  2004/08/06 12:37:25  roboos
% fixed minor bugs
%
% Revision 1.1  2004/08/06 08:55:23  roboos
% new implementation
%

% these are suitable RGB colors
skin   = [255 213 119]/255;
brain  = [202 100 100]/255;
cortex = [255 213 119]/255;

% set the defaults
if ~isfield(cfg, 'surface_facecolor'), cfg.surface_facecolor = skin;   end
if ~isfield(cfg, 'surface_edgecolor'), cfg.surface_edgecolor = 'none'; end
if ~isfield(cfg, 'surface_facealpha'), cfg.surface_facealpha = 0.7;    end

if ~isfield(cfg, 'vol') && ~isfield(cfg, 'hdmfile')
  cfg.vol = [];  % FIXME why is this empty setting neccessary?
end

if nargin<2
  data = [];
end

% set the defaults that apply both to EEG and MEG
if ~isfield(cfg, 'spheremesh'),       cfg.spheremesh = 642;          end
if ~isfield(cfg, 'plotsensors'),      cfg.plotsensors = 'yes';       end
if ~isfield(cfg, 'plotheadsurface'),  cfg.plotheadsurface = 'yes';   end
if ~isfield(cfg, 'plotgrid'),         cfg.plotgrid = 'yes';          end
if ~isfield(cfg, 'plotinside'),       cfg.plotinside = 'yes';        end
if ~isfield(cfg, 'plotoutside'),      cfg.plotoutside = 'no';        end
if ~isfield(cfg, 'plotheadshape'),    cfg.plotheadshape = 'yes';     end

% extract/read the gradiometer and volume conductor
[vol, sens, cfg] = prepare_headmodel(cfg, data);

% extract/construct the dipole grid
if strcmp(cfg.plotgrid, 'yes') && isfield(cfg, 'grid')
  sourcegrid = prepare_dipole_grid(cfg, vol, sens);
else
  sourcegrid     = [];
  sourcegrid.pos = zeros(0,3);
  sourcegrid.inside = [];
  sourcegrid.outside = [];
end

% determine the type of input data
ismeg          = isfield(sens, 'ori');
iseeg          = (~ismeg);
isbem          = isfield(vol, 'bnd');
issphere       = isfield(vol, 'r');
ismultisphere  = isfield(vol, 'r') && length(vol.r)>4;
issinglesphere = isfield(vol, 'r') && length(vol.r)==1;

if ismeg
  % sensors describe MEG data, set the corresponding defaults
  if ~isfield(cfg, 'plotspheres'),      cfg.plotspheres = 'no';        end
  if ~isfield(cfg, 'plotspherecenter'), cfg.plotspherecenter = 'no';   end
  if ~isfield(cfg, 'plotlines'),        cfg.plotlines = 'yes';         end
  if ~isfield(cfg, 'plotcoil'),         cfg.plotcoil = 'no';           end
  if ~isfield(cfg, 'headshape'),        cfg.headshape = [];            end
else
  % sensors describe EEG data, set the corresponding defaults
  if ~isfield(cfg, 'plotspheres'),      cfg.plotspheres = 'no';        end
  if ~isfield(cfg, 'plotspherecenter'), cfg.plotspherecenter = 'no';   end
  if ~isfield(cfg, 'plotlines'),        cfg.plotlines = 'yes';         end
end

if     ismeg && strcmp(sensortype(sens), 'ctf151')
  Nsensors = length(channelselection('MEG', sens.label));
elseif ismeg && strcmp(sensortype(sens), 'ctf275')
  Nsensors = length(channelselection('MEG', sens.label));
else
  % assume that the number of sensors to be plotted is the same as the number of channels
  % this implies that channel 1:Nsensors corresponds with the bottom coils
  % of a gradiometer MEG system (which is the case for a CTF system by
  % virtue of ctf2grad).
  Nsensors = length(sens.label);
end

if issphere
  % determine the number of spheres in the volume model
  Nspheres = length(vol.r);
  fprintf('spherical head model with %d spheres\n', Nspheres);
end

clf
hold on
axis equal
axis vis3d
if ismeg
  axis off
else
  axis on
  grid on
  xlabel('x')
  ylabel('y')
  zlabel('z')
end

if iseeg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plotting for EEG
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if strcmp(cfg.plotgrid, 'yes')
    if strcmp(cfg.plotinside, 'yes')
      plot3(sourcegrid.pos(sourcegrid.inside,1), sourcegrid.pos(sourcegrid.inside,2), sourcegrid.pos(sourcegrid.inside,3), 'k.');
    end
    if strcmp(cfg.plotoutside, 'yes')
      plot3(sourcegrid.pos(sourcegrid.outside,1), sourcegrid.pos(sourcegrid.outside,2), sourcegrid.pos(sourcegrid.outside,3), 'k.');
    end
  end % plotgrid

  if strcmp(cfg.plotsensors, 'yes')
    plot3(sens.pnt(1:Nsensors,1), sens.pnt(1:Nsensors,2), sens.pnt(1:Nsensors,3), 'b*');
  end % plotsensors

  if strcmp(cfg.plotheadsurface, 'yes')  && ~isempty(vol)
    [pnt, tri] = headsurface(vol, sens);
    h = patch('Faces', tri, 'Vertices', pnt);
    set(h, 'facecolor', cfg.surface_facecolor);
    set(h, 'edgecolor', cfg.surface_edgecolor);
    set(h, 'facealpha', cfg.surface_facealpha);
  end % plotheadsurface

  if strcmp(cfg.plotlines, 'yes') && ~isempty(vol)
    if isbem
      % project the electrodes on the skin surface, on the nearest triangle
      [el] = project_elec(sens.pnt, vol.bnd(vol.skin).pnt, vol.bnd(vol.skin).tri);
      % this returns [tri, la, mu]
      tri = el(:,1);
      la  = el(:,2);
      mu  = el(:,3);
      for i=1:Nsensors
        v1  = vol.bnd(vol.skin).pnt(vol.bnd(vol.skin).tri(tri(i),1),:);
        v2  = vol.bnd(vol.skin).pnt(vol.bnd(vol.skin).tri(tri(i),2),:);
        v3  = vol.bnd(vol.skin).pnt(vol.bnd(vol.skin).tri(tri(i),3),:);
        prj(i,:) = routlm(v1, v2, v3, la(i), mu(i));
      end
    elseif issphere
      % project the electrodes onto the sphere surface, towards the origin
      nrm = sqrt(sum(sens.pnt.^2,2));
      prj = vol.r(vol.skin) * sens.pnt ./ [nrm nrm nrm];
    else
      % in case that no known volume conductor is specified
      prj = sens.pnt;
    end
    x = [sens.pnt(:,1) prj(:,1)]';
    y = [sens.pnt(:,2) prj(:,2)]';
    z = [sens.pnt(:,3) prj(:,3)]';
    line(x, y, z, 'color', 'm');
  end % plotlines

elseif ismeg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plotting for MEG
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if strcmp(cfg.plotgrid, 'yes')
    if strcmp(cfg.plotinside, 'yes')
      plot3(sourcegrid.pos(sourcegrid.inside,1), sourcegrid.pos(sourcegrid.inside,2), sourcegrid.pos(sourcegrid.inside,3), 'k.');
    end
    if strcmp(cfg.plotoutside, 'yes')
      plot3(sourcegrid.pos(sourcegrid.outside,1), sourcegrid.pos(sourcegrid.outside,2), sourcegrid.pos(sourcegrid.outside,3), 'k.');
    end
  end % plotgrid

  if strcmp(cfg.plotsensors, 'yes')
    % plot only the bottom coil
    plot3(sens.pnt(1:Nsensors,1), sens.pnt(1:Nsensors,2), sens.pnt(1:Nsensors,3), 'b*');
  end % plotsensors

  if strcmp(cfg.plotcoil, 'yes')
    pnt = sens.pnt;
    ori = sens.ori;
    x = [pnt(:,1) pnt(:,1)+ori(:,1)]';
    y = [pnt(:,2) pnt(:,2)+ori(:,2)]';
    z = [pnt(:,3) pnt(:,3)+ori(:,3)]';
    line(x, y, z, 'color', 'g');
    plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'g.');
  end % plotsensors

  if ~isempty(cfg.headshape) && strcmp(cfg.plotheadshape, 'yes')
    % get the surface describing the head shape
    if isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
      % use the headshape surface specified in the configuration
      headshape.pnt = cfg.headshape.pnt;
    elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
      % use the headshape points specified in the configuration
      headshape.pnt = cfg.headshape;
    elseif ischar(cfg.headshape) && filetype(cfg.headshape, 'ctf_shape')
      % read the headshape from file
      headshape = read_ctf_shape(cfg.headshape);
    elseif ischar(cfg.headshape) && filetype(cfg.headshape, '4d_hs')
      % read the headshape from file
      headshape     = []; 
      headshape.pnt = read_bti_hs(cfg.headshape);
    else
      error('cfg.headshape is not specified correctly')
    end
    plot3(headshape.pnt(:,1), headshape.pnt(:,2), headshape.pnt(:,3), 'r.');
  end % plotheadshape

  if strcmp(cfg.plotspheres, 'yes') && ~isempty(vol)
    % create a triangulated unit sphere
    if cfg.spheremesh==42
      [pnt0, tri] = icosahedron42;
    elseif cfg.spheremesh==162
      [pnt0, tri] = icosahedron162;
    elseif cfg.spheremesh==642
      [pnt0, tri] = icosahedron642;
    end
    Nvertices = size(pnt0,1);
    for i=1:Nspheres
      % scale and shift the unit sphere to the proper location
      pnt = pnt0*vol.r(i) + repmat(vol.o(i,:),Nvertices,1);
      if Nspheres>4
        distance = sqrt(sum((pnt - repmat(sens.pnt(i,:),Nvertices,1)).^2, 2));
        distance = (distance-min(distance))/range(distance);
      else
        distance = zeros(Nvertices, 1);
      end
      h = patch('Faces', tri, 'Vertices', pnt);
      % set(h, 'FaceVertexCData', 0.5*ones(length(distance),30));
      % set(h, 'FaceVertexCData', [0 0 1]);
      set(h, 'FaceVertexAlphaData', (1-distance).*(distance<0.1));
      % set(h, 'AlphaDataMapping', 'direct');
      set(h, 'edgealpha', 'interp')
      % set(h, 'edgealpha', 1)
      set(h, 'facealpha', 'interp')
      % set(h, 'facealpha', 0)
      set(h, 'facecolor', 'none');
      % set(h, 'linestyle', 'none');
    end
  end % plotspheres

  if strcmp(cfg.plotspherecenter, 'yes') && ~isempty(vol)
    plot3(vol.o(:,1), vol.o(:,2), vol.o(:,3), 'k.');
  end % plotspherecenter

  if strcmp(cfg.plotheadsurface, 'yes') && ~isempty(vol)
    % estimate the head surface from the spheres and gradiometers
    [pnt, tri] = headsurface(vol, sens);
    h = patch('Faces', tri, 'Vertices', pnt);
    set(h, 'facecolor', cfg.surface_facecolor);
    set(h, 'edgecolor', cfg.surface_edgecolor);
    set(h, 'facealpha', cfg.surface_facealpha);
  end % plotheadsurface

  if strcmp(cfg.plotlines, 'yes') && ismultisphere
    % first determine the indices of the relevant gradiometers.
    % [sel_g, sel_v] = match_str(sens.label, vol.label);
    sel_g = 1:Nsensors;
    sel_v = 1:Nsensors;
    % create head-surface points from multisphere head-model.
    dir     = sens.pnt(sel_g,:) - vol.o(sel_v,:);
    dist    = sqrt(sum(dir.*dir,2));
    pnt0    = repmat((vol.r(sel_v)./dist),1,3).*dir + vol.o(sel_v,:);
    pnt1    = sens.pnt(sel_g,:);
    f = gcf;
    figure
    hist(sqrt(sum((pnt0-pnt1).^2, 2)));
    xlabel('distance (a.u.)');
    title('histogram of distances from the bottom coil to the head surface');
    figure(f);
    x = [pnt0(:,1) pnt1(:,1)]';
    y = [pnt0(:,2) pnt1(:,2)]';
    z = [pnt0(:,3) pnt1(:,3)]';
    line(x, y, z, 'color', 'm');
  end % plotlines

end % iseeg or ismeg

lighting gouraud
camlight
hold off

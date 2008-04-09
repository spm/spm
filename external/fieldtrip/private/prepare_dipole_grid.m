function [grid, cfg] = prepare_dipole_grid(cfg, vol, sens);

% PREPARE_DIPOLE_GRID helps to make a regular grid with dipoles that can be
% used for scanning, e.g. as initial search for dipole fitting or for
% beamforming.
%
% Configuration options for generating a regular 3-D grid
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
%
% Configuration options for a predefined grid
%   cfg.grid.pos        = N*3 matrix with position of each source
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%   cfg.grid.inside     = vector with indices of the sources inside the brain (optional)
%   cfg.grid.outside    = vector with indices of the sources outside the brain (optional)
% The following fields are not used in this function, but will be copied along to the output
%   cfg.grid.leadfield
%   cfg.grid.filter or alternatively cfg.grid.avg.filter
%   cfg.grid.subspace
%   cfg.grid.lbex
%
% Configuration options for cortex segmentation, i.e. for placing dipoles in grey matter
%   cfg.mri           = can be filename, MRI structure or segmented MRI structure
%   cfg.sourceunits   = 'mm' or 'cm' (default is 'cm')
%   cfg.mriunits      = 'mm' or 'cm' (default is 'mm')
%   cfg.threshold     = 0.1, relative to the maximum value in the segmentation
%   cfg.smooth        = 5, smoothing in voxels
%
% Other configuration options
%   cfg.headshape   = filename of headshape (optional)
%   cfg.inwardshift = depth of the bounding layer for the source space, relative to the head model surface (default = 0)
%   cfg.symmetry    = 'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
%   cfg.tightgrid   = 'yes' or 'no' (default is automatic)
%
% See also SOURCEANALYSIS, DIPOLEFITTING, PREPARE_LEADFIELD

% Copyright (C) 2004-2007, Robert Oostenveld
%
% $Log: prepare_dipole_grid.m,v $
% Revision 1.35  2007/12/11 11:12:44  roboos
% remember dipole moment if specified
%
% Revision 1.34  2007/06/14 10:25:51  jansch
% changed order of handling of haspos and hasgrid in order to be able to create
% mni-based grids
%
% Revision 1.33  2007/06/07 13:00:14  roboos
% fixed bug in detection of hasgrid, thanks to Juan
%
% Revision 1.32  2007/05/16 12:22:00  roboos
% fixed some small bugs
%
% Revision 1.31  2007/05/16 11:29:11  roboos
% included the generation of surface grids based on headshape or volume (used in megrealign)
%
% Revision 1.30  2007/01/03 15:04:44  roboos
% default cfg.tightgrid=no if cfg.grid.pos specified
%
% Revision 1.29  2007/01/02 13:35:16  roboos
% added option cfg.tightgrid to make the tightening accessible from the outside
%
% Revision 1.28  2006/10/12 09:58:15  jansch
% added default for cfg.symmetry
%
% Revision 1.27  2006/10/12 08:44:48  roboos
% implemented construction of symmetric dipole pair, cfg.symmetry
%
% Revision 1.26  2006/07/24 08:22:20  roboos
% added option cfg.inwardshift, default=0 (can be negative for an outward shift)
%
% Revision 1.25  2006/07/05 10:21:10  roboos
% moved cfg.xgrid/ygrid/zgrid/resolution options into cfg.grid substructure
% added code for automatic estimation of source units in case cfg.grid.resolution is unspecified
% cleaned up code a little bit, cleaned up documentation
%
% Revision 1.24  2006/06/08 07:50:57  roboos
% updated the conversion between the source and MRI units (support mm,cm,dm,m for both)
%
% Revision 1.23  2006/05/23 10:16:11  roboos
% Also update the cfg in case of tightening the grid.

% set the defaults
if ~isfield(cfg, 'spheremesh'),       cfg.spheremesh = 642;       end
if ~isfield(cfg, 'symmetry'),         cfg.symmetry = [];          end
if ~isfield(cfg, 'grid'),             cfg.grid = [];              end

if isfield(cfg, 'grid') && isempty(cfg.grid)
  % for backward compatibility with old cfg, move the options to substructure
  if isfield(cfg, 'resolution'),      cfg.grid.resolution = cfg.resolution; cfg = rmfield(cfg, 'resolution'); end
  if isfield(cfg, 'xgrid'),           cfg.grid.xgrid = cfg.xgrid; cfg = rmfield(cfg, 'xgrid'); end
  if isfield(cfg, 'ygrid'),           cfg.grid.ygrid = cfg.ygrid; cfg = rmfield(cfg, 'ygrid'); end
  if isfield(cfg, 'zgrid'),           cfg.grid.zgrid = cfg.zgrid; cfg = rmfield(cfg, 'zgrid'); end
end

if isfield(cfg, 'grid') && ~isempty(cfg.grid)
  % generate an automatic grid that fits around the sensor array
  if ~isfield(cfg.grid, 'xgrid'),     cfg.grid.xgrid = 'auto';    end
  if ~isfield(cfg.grid, 'ygrid'),     cfg.grid.ygrid = 'auto';    end
  if ~isfield(cfg.grid, 'zgrid'),     cfg.grid.zgrid = 'auto';    end
  if ~isfield(cfg, 'inwardshift'),    cfg.inwardshift = 0;        end % in this case for inside detection
end

if isfield(cfg, 'mri')
  % set the defaults for cortex segmentation
  if ~isfield(cfg, 'threshold'),      cfg.threshold = 0.1;        end % relative
  if ~isfield(cfg, 'smooth');         cfg.smooth    = 5;          end % in voxels
  if ~isfield(cfg, 'sourceunits');    cfg.sourceunits = 'cm';     end
  if ~isfield(cfg, 'mriunits');       cfg.mriunits = 'mm';        end
end

% there are a number of ways in which the grid can be constructed
hasgrid  = isfield(cfg.grid, 'resolution') || isfield(cfg.grid, 'xgrid');  % from the specified configuration (xgrid/ygrid/zgrid/resolution)
haspos   = isfield(cfg.grid, 'pos');                                       % the dipole positions are fully specified
hasmom   = isfield(cfg.grid, 'mom');                                       % the dipole moments/orientations are fully specified
hasmri   = isfield(cfg, 'mri');                                            % a segmented cortex will be used to determine positions
hasshape = isfield(cfg, 'headshape') && ~isempty(cfg.headshape);           % construct grid from the inward shifted headshape

if ~isfield(cfg, 'tightgrid')
  % If any of the grid options is 'auto, an initial grid is constructed
  % that spans a box which completely encompasses the sensor array.
  % Such a grid is much larger than what is needed to span the brain
  % compartment, therefore it will be tightened around the brain
  % compartment at the end.
  if hasgrid
    if isfield(cfg.grid, 'pos')
      cfg.tightgrid = 'no';
    elseif ischar(cfg.grid.xgrid) && strcmp(cfg.grid.xgrid, 'auto')
      cfg.tightgrid = 'yes';
    elseif ischar(cfg.grid.ygrid) && strcmp(cfg.grid.ygrid, 'auto')
      cfg.tightgrid = 'yes';
    elseif ischar(cfg.grid.zgrid) && strcmp(cfg.grid.zgrid, 'auto')
      cfg.tightgrid = 'yes';
    else
      cfg.tightgrid = 'no';
    end
  else
    cfg.tightgrid = 'no';
  end
end

% start with an empty grid
grid = [];

if hasmri
  % convert the source/functional data into the same units as the anatomical MRI
  scale = 1;
  switch cfg.sourceunits
    case 'mm'
      scale = scale / 1000;
    case 'cm'
      scale = scale / 100;
    case 'dm'
      scale = scale / 10;
    case 'm'
      scale = scale / 1;
    otherwise
      error('unknown physical dimension in cfg.sourceunits');
  end
  switch cfg.mriunits
    case 'mm'
      scale = scale * 1000;
    case 'cm'
      scale = scale * 100;
    case 'dm'
      scale = scale * 10;
    case 'm'
      scale = scale * 1;
    otherwise
      error('unknown physical dimension in cfg.mriunits');
  end
end

if haspos
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % a grid is already specified in the configuration, reuse as much of the
  % prespecified grid as possible (but only known objects)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  grid.pos = cfg.grid.pos;
  if hasmom
    grid.mom = cfg.grid.mom;
  end
  if isfield(cfg.grid, 'xgrid')
    grid.xgrid = cfg.grid.xgrid;
  end
  if isfield(cfg.grid, 'ygrid')
    grid.ygrid = cfg.grid.ygrid;
  end
  if isfield(cfg.grid, 'zgrid')
    grid.zgrid = cfg.grid.zgrid;
  end
  if isfield(cfg.grid, 'dim')
    grid.dim = cfg.grid.dim;
  elseif isfield(grid, 'xgrid') && isfield(grid, 'ygrid') && isfield(grid, 'zgrid')
    grid.dim = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
  end
  if isfield(cfg.grid, 'inside')
    grid.inside = cfg.grid.inside;
  end
  if isfield(cfg.grid, 'outside')
    grid.outside = cfg.grid.outside;
  end
  if isfield(cfg.grid, 'lbex')
    grid.lbex = cfg.grid.lbex;
  end
  if isfield(cfg.grid, 'subspace')
    grid.subspace = cfg.grid.subspace;
  end
  if isfield(cfg.grid, 'leadfield')
    grid.leadfield = cfg.grid.leadfield;
  end
  if isfield(cfg.grid, 'filter')
    grid.filter = cfg.grid.filter;
  elseif isfield(cfg.grid, 'avg')
    % the precomputed filters could also be in the average sub-structure
    if isfield(cfg.grid.avg, 'filter')
      grid.filter = cfg.grid.avg.filter;
    end
  end
  % this will remove all remaining unknown objects from the configuration
  cfg.grid = grid;

elseif hasgrid
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % construct a regular 3D grid that spans a box encompassing all electrode
  % or gradiometer coils, this will typically also cover the complete brain
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (strcmp(cfg.grid.xgrid, 'auto') || strcmp(cfg.grid.ygrid, 'auto') || strcmp(cfg.grid.zgrid, 'auto')) && ~isfield(cfg.grid, 'resolution')
    % the units of the head coordinate system are unknown and can be mm, cm or m
    % try to determine the units by looking at the spatial extent of the sensor array
    switch floor(log10(max(range(sens.pnt))))
      case -1
        % cfg.sourceunits = 'm';
        cfg.grid.resolution = 0.01;
      case 0
        % cfg.sourceunits = 'dm';
        cfg.grid.resolution = 0.1;
      case 1
        % cfg.sourceunits = 'cm';
        cfg.grid.resolution = 1;
      case 2
        % cfg.sourceunits = 'mm';
        cfg.grid.resolution = 10;
    end
    fprintf('using an automatic grid resolution of %d\n', cfg.grid.resolution);
  end
  if ischar(cfg.grid.xgrid) && strcmp(cfg.grid.xgrid, 'auto')
    cfg.grid.xgrid = floor(min(sens.pnt(:,1))):cfg.grid.resolution:ceil(max(sens.pnt(:,1)));
  end
  if ischar(cfg.grid.ygrid) && strcmp(cfg.grid.ygrid, 'auto')
    cfg.grid.ygrid = floor(min(sens.pnt(:,2))):cfg.grid.resolution:ceil(max(sens.pnt(:,2)));
  end
  if ischar(cfg.grid.zgrid) && strcmp(cfg.grid.zgrid, 'auto')
    cfg.grid.zgrid = floor(min(sens.pnt(:,3))):cfg.grid.resolution:ceil(max(sens.pnt(:,3)));
  end
  % create the regular 3-D dipole grid
  grid.xgrid = cfg.grid.xgrid;
  grid.ygrid = cfg.grid.ygrid;
  grid.zgrid = cfg.grid.zgrid;
  grid.dim   = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
  [X, Y, Z]  = ndgrid(cfg.grid.xgrid, cfg.grid.ygrid, cfg.grid.zgrid);
  grid.pos   = [X(:) Y(:) Z(:)];

elseif hasmri
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % construct a grid based on the segmented MRI that is provided in the
  % configuration, only voxels in gray matter will be used
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isfield(cfg.grid, 'resolution')
    switch cfg.sourceunits
      case 'mm'
        cfg.grid.resolution = 10;
      case 'cm'
        cfg.grid.resolution = 1;
      case 'dm'
        cfg.grid.resolution = 0.1;
      case 'm'
        cfg.grid.resolution = 0.01;
    end
  end
  if isstr(cfg.mri)
    % read the segmentation from file
    mri      = read_fcdc_mri(cfg.mri);
    mri.gray = double(mri.anatomy);
  elseif isstruct(cfg.mri) && ~isfield(cfg.mri, 'gray')
    % looks like a segmentation that has already been loaded in memory
    mri      = cfg.mri;
    mri.gray = double(mri.anatomy);
  elseif isstruct(cfg.mri) && isfield(cfg.mri, 'gray')
    % looks like a complete segmentation from VOLUMESEGMENT
    mri      = cfg.mri;
    mri.gray = double(mri.gray);
  else
    error('cannot determine the format of the segmentation in cfg.mri');
  end

  % for compatibility with old MRI volumes
  mri = fixvolume(mri);

  % apply a smoothing of a certain amount of voxels
  if ~isstr(cfg.smooth) && cfg.smooth>1
    fprintf('smoothing gray matter segmentation with %d voxels\n', cfg.smooth);
    mri.gray = spm_conv(mri.gray, cfg.smooth);
  end

  % determine for each voxel whether it belongs to the cortex
  if isfield(cfg, 'threshold')
    fprintf('thresholding gray matter segmentation at a relative value of %f\n', cfg.threshold);
    head = mri.gray./max(mri.gray(:)) > cfg.threshold;
  else
    error('you must specify cfg.threshold for cortex segmentation');
  end

  ind                 = find(head(:));
  fprintf('%d from %d voxels in the segmentation are marked as cortex (%.0f%%)\n', length(ind), prod(size(head)), 100*length(ind)/prod(size(head)));
  [X,Y,Z]             = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));   % create the grid in MRI-coordinates
  posmri              = [X(ind) Y(ind) Z(ind) ones(length(ind),1)];         % take only the inside voxels
  poshead             = mri.transform * posmri';                            % transform to head coordinates
  poshead             = poshead(1:3,:)';
  posmri              = posmri(:,1:3);
  resolution          = cfg.grid.resolution*scale;                                        % source and mri can be expressed in different units (e.g. cm and mm)
  xgrid               = floor(min(poshead(:,1))):resolution:ceil(max(poshead(:,1)));      % create the grid in head-coordinates
  ygrid               = floor(min(poshead(:,2))):resolution:ceil(max(poshead(:,2)));      % with 'consistent' x,y,z definitions
  zgrid               = floor(min(poshead(:,3))):resolution:ceil(max(poshead(:,3)));
  [X,Y,Z]             = ndgrid(xgrid,ygrid,zgrid);
  pos2head            = [X(:) Y(:) Z(:) ones(length(X(:)),1)]';
  pos2mri             = mri.transform \ pos2head;                                         % transform to MRI-coordinates
  pos2mri             = round(pos2mri(1:3,:))';
  pos2head            = pos2head(1:3,:)';
  pos2mri             = pos2mri(:,1:3);
  % it might be that the box with the points does not completely fit into the MRI
  sel = find(pos2mri(:,1)<1 |  pos2mri(:,1)>size(head,1) | ...
    pos2mri(:,2)<1 |  pos2mri(:,2)>size(head,2) | ...
    pos2mri(:,3)<1 |  pos2mri(:,3)>size(head,3));
  if isempty(sel)
    % use the efficient implementation
    inside = head(sub2ind(mri.dim, pos2mri(:,1), pos2mri(:,2), pos2mri(:,3)));
  else
    % only loop over the points that can be dealt with
    inside = zeros(length(xgrid)*length(ygrid)*length(zgrid), 1);
    for i=setdiff(1:size(pos2mri,1), sel(:)')
      inside(i) = head(pos2mri(i,1), pos2mri(i,2), pos2mri(i,3));
    end
  end
  inside = find(inside);

  grid.pos            = pos2head/scale;                                     % convert to source units
  grid.xgrid          = xgrid/scale;                                        % convert to source units
  grid.ygrid          = ygrid/scale;                                        % convert to source units
  grid.zgrid          = zgrid/scale;                                        % convert to source units
  grid.dim            = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
  grid.inside         = inside(:);
  grid.outside        = setdiff(1:size(grid.pos,1),grid.inside)';

  fprintf('the regular 3D grid encompassing the cortex contains %d grid points\n', size(grid.pos,1));
  fprintf('%d grid points inside gray matter\n', length(grid.inside));
  fprintf('%d grid points outside gray matter\n', length(grid.outside));

elseif hasshape
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the headshape  to make a superficial dipole layer (e.g.
  % for megrealign). Assume that all points are inside the volume.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  grid.pos = headsurface([], [], 'headshape', cfg.headshape, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  grid.inside  = 1:size(grid.pos,1);
  grid.outside = [];

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the volume conduction model to make a superficial dipole layer (e.g.
  % for megrealign). Assume that all points are inside the volume.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  grid.pos = headsurface(vol, sens, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  grid.inside  = 1:size(grid.pos,1);
  grid.outside = [];

end

% determine the dipole locations inside the brain volume
if ~isfield(grid, 'inside') && ~isfield(grid, 'outside')
  if isempty(vol)
    % an empty vol in combination with gradiometers indicates a magnetic dipole
    % in an infinite vacuum, i.e. all dipoles can be considered to be inside
    grid.inside = 1:size(grid.pos,1);
    grid.outside = [];
  else
    if isfield(sens, 'ori') && isfield(sens, 'pnt') && isfield(sens, 'tra')
      % in case of MEG, make a triangulation of the outermost surface
      if isfield(cfg, 'headshape')
        % use the specified headshape to construct the bounding triangulation
        [pnt, tri] = headsurface(vol, sens, 'headshape', cfg.headshape, 'inwardshift', cfg.inwardshift, 'surface', 'skin');
      else
        % use the volume conductor model to construct the bounding triangulation
        [pnt, tri] = headsurface(vol, sens, 'inwardshift', cfg.inwardshift, 'surface', 'skin');
      end
    else
      % in case of EEG, make a triangulation of the innermost surface
      [pnt, tri] = headsurface(vol, sens, 'inwardshift', cfg.inwardshift, 'surface', 'brain');
    end
    % determine the dipole positions that are inside the triangulated surface
    tmp = bounding_mesh(grid.pos, pnt, tri);
    grid.inside  = find(tmp==1);
    grid.outside = find(tmp==0);
  end
elseif ~isfield(grid, 'inside')
  grid.inside = setdiff(1:size(grid.pos,1), grid.outside);
elseif ~isfield(grid, 'outside')
  grid.outside = setdiff(1:size(grid.pos,1), grid.inside);
end

if strcmp(cfg.tightgrid, 'yes')
  fprintf('%d dipoles inside, %d dipoles outside brain\n', length(grid.inside), length(grid.outside));
  fprintf('making tight grid\n');
  xmin = min(grid.pos(grid.inside,1));
  ymin = min(grid.pos(grid.inside,2));
  zmin = min(grid.pos(grid.inside,3));
  xmax = max(grid.pos(grid.inside,1));
  ymax = max(grid.pos(grid.inside,2));
  zmax = max(grid.pos(grid.inside,3));
  xmin_indx = find(grid.xgrid==xmin);
  ymin_indx = find(grid.ygrid==ymin);
  zmin_indx = find(grid.zgrid==zmin);
  xmax_indx = find(grid.xgrid==xmax);
  ymax_indx = find(grid.ygrid==ymax);
  zmax_indx = find(grid.zgrid==zmax);
  sel =       (grid.pos(:,1)>=xmin & grid.pos(:,1)<=xmax); % select all grid positions inside the tight box
  sel = sel & (grid.pos(:,2)>=ymin & grid.pos(:,2)<=ymax); % select all grid positions inside the tight box
  sel = sel & (grid.pos(:,3)>=zmin & grid.pos(:,3)<=zmax); % select all grid positions inside the tight box
  grid.pos   = grid.pos(sel,:);
  % update the inside and outside vector
  tmp = zeros(1,prod(grid.dim));
  tmp(grid.inside)  = 1;        % these are originally inside the brain
  tmp(grid.outside) = 0;        % these are originally outside the brain
  tmp               = tmp(sel); % within the tight box, these are inside the brain
  grid.inside  = find(tmp);
  grid.outside = find(~tmp);
  grid.xgrid   = grid.xgrid(xmin_indx:xmax_indx);
  grid.ygrid   = grid.ygrid(ymin_indx:ymax_indx);
  grid.zgrid   = grid.zgrid(zmin_indx:zmax_indx);
  grid.dim     = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
end

fprintf('%d dipoles inside, %d dipoles outside brain\n', length(grid.inside), length(grid.outside));

% apply the symmetry constraint, i.e. add a symmetric dipole for each location defined sofar
% set up the symmetry constraints
if ~isempty(cfg.symmetry)
  if strcmp(cfg.symmetry, 'x')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 -1 1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 -x1 y
  elseif strcmp(cfg.symmetry, 'y')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 1 -1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 -y
  elseif strcmp(cfg.symmetry, 'z')
    reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    mirror = [1 1 1 1 1 -1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 y1
  else
    error('unrecognized symmetry constraint');
  end
  fprintf('each source describes two dipoles with symmetry along %s axis\n', cfg.symmetry);
  % expand the number of parameters from one (3) to two dipoles (6)
  grid.pos = grid.pos(:,expand) .* repmat(mirror, size(grid.pos,1), 1);
end

% update the configuration for consistency
cfg.grid = grid;

function sourceplot(cfg, interp)

% SOURCEPLOT displays three orthogonal slices of anatomical and/or functional
% volume data. If an anatomical volume is present in the input, it will
% always be plotted. The source parameter that is specified in the
% configuration will be plotted as an overlay on the anatomy, where the
% mask parameter can be used to determine the opacity.
%
% Use as
%   sourceplot(cfg, functional)
% where the functional data is the output of SOURCEANALYSIS,
% SOURCESTATISTICS, SOURCEINTERPOLATE or NORMALISEVOLUME and 
% the configuration contains any of the following fields
%   cfg.location            = 'min', 'max', 'interactive' or [x y z]
%   cfg.locationcoordinates = 'head' or 'voxel' (default = 'head')
%   cfg.funparameter        = string with the functional parameter of interest
%   cfg.maskparameter       = string with an optional mask parameter
%   cfg.maskcolmin          = mask value mapped to the lowest opacity, i.e. completely transparent (default ='auto')
%   cfg.maskcolmax          = mask value mapped to the highest opacity, i.e. non-transparent (default = 'auto')
%   cfg.colmin              = source value mapped to the lowest color (default = 'auto')
%   cfg.colmax              = source value mapped to the highest color (default = 'auto')
%
% See also SOURCEANALYSIS, SOURCESTATISTICS, SOURCEINTERPOLATE

% Undocumented local options:
% cfg.TTlookup      = 'yes' or 'no' (default)
% cfg.TTqueryrange  = number (default 3)
% cfg.colormap
% cfg.axis          = 'on' (default) or 'off'
% cfg.crosshair     = 'on' (default) or 'off'
% cfg.colorbar      = 'on' (default) or 'off'
% cfg.location      = 'maxproject', creates a glass-brain representation with an outline of a template brain only works with 1mm resolution normalised images! ask JM

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: sourceplot_old.m,v $
% Revision 1.1  2007/02/08 12:34:26  roboos
% this is a copy of the sourcelot function just prior to Ingrid rewriting it
%
% Revision 1.38  2006/09/18 14:04:45  jansch
% implemented maxproject as an undocumented option for cfg.location. This gives
% a 'glass-brain' projection of the functional data
%
% Revision 1.37  2006/07/13 08:48:46  ingnie
% fixed typo's in documentation
%
% Revision 1.36  2006/05/31 07:07:15  roboos
% updated documentation
%
% Revision 1.35  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.34  2006/03/29 16:14:44  roboos
% do a cla (clear axes) in each subplot when plotting interactive (thanks to Vladimir)
%
% Revision 1.33  2006/03/07 08:33:29  roboos
% removed debugging keyboard statement
%
% Revision 1.32  2006/03/07 08:32:42  roboos
% removed mm units from the print statement about the location on screen
% removed the cfg.flipdim option, it invalidated the reported position
%
% Revision 1.31  2006/03/02 13:21:25  jansch
% downgraded version 1.30 to version 1.27 since the newer version was buggy
%
% Revision 1.27  2006/01/30 17:11:04  roboos
% added cfg.flipdim option, flips the 3rd dimension
%
% Revision 1.26  2006/01/25 11:01:04  roboos
% add xgrid/ygrid/zgrid in voxel coordinates
%
% Revision 1.25  2006/01/24 21:26:09  roboos
% fixed tmpcfg.parameter for downsample
% moved parameterselection up in the code
%
% Revision 1.24  2006/01/24 14:32:23  roboos
% list cfg.TTqueryrange and cfg.TTlookup as undocumented options (not in main help, but just below it)
%
% Revision 1.23  2006/01/05 13:35:41  roboos
% Changed parameterselection subfunction, select the first element
% from the cell-array.
% Do not add xgrid/ygrid/zgrid or transform, call the private function
% grid2transform to ensure that the volume is described using the
% homogenous coordinate transformation matrix.
% Always pass the volume to the downsamplevolume function, that
% function is smart enough to notice cfg.downsample=1.
%
% Revision 1.22  2005/10/14 15:48:01  roboos
% added support for cfg.maskcolmin and cfg.maskcolmax, similar to sliceinterp
%
% Revision 1.21  2005/08/19 16:57:48  roboos
% changed the default opacity to 0.5
% use the new subfunction parameterselection() for selecting the functional and mask volume
%
% Revision 1.20  2005/08/17 19:22:09  roboos
% removed obsolete reference to cfg.TTmask
%
% Revision 1.19  2005/07/29 13:23:10  roboos
% renamed searchrange to queryrange
% added coordinate translation using mni2tal()
%
% Revision 1.18  2005/07/20 15:32:28  roboos
% added experimental code to print the anatomical labels from the TT atlas
%
% Revision 1.17  2005/06/17 10:59:38  roboos
% add x/y/zgrid if not present in input data
%
% Revision 1.16  2005/05/17 17:50:38  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.15  2005/03/07 17:09:21  roboos
% added a custom colorbar for the functional data in the 4th subplot (thanks to Markus)
%
% Revision 1.14  2005/03/07 15:07:59  roboos
% added option cfg.locationcoordinates, can be 'head' (default) or 'voxel'
%
% Revision 1.13  2005/03/03 10:38:26  roboos
% changed axes labels from x/y/z into i/j/k, since axes correspond with voxel indices
%
% Revision 1.12  2005/02/28 17:58:59  roboos
% fixed bug in selection of empty funparameter
%
% Revision 1.11  2005/02/21 12:29:08  roboos
% removed incorrect inv() from interp.transform when interactive plotting
%
% Revision 1.10  2005/02/11 17:59:01  roboos
% fixed bug in computation of headcoordinates (removed inv)
%
% Revision 1.9  2005/02/11 14:49:11  roboos
% moved clipping and scaling from subfunction to mainfunction
% implemented cfg.downsample
%
% Revision 1.8  2005/02/11 13:28:53  roboos
% completely new implementation, basically based upon misc/volplot and on Nienke's singleplot
%
% Revision 1.7  2004/08/26 10:48:25  roboos
% modified to make it consistent again with the new volume (source+mri) structure
% removed dependency on warp3d, implemented homogenous transform directly

% sometimes it is desirable to specify the cfg as a cell array
if iscell(cfg)
  var = cfg(1:2:length(cfg)); % select the odd elements
  val = cfg(2:2:length(cfg)); % select the even elements
  % convert the configuration to a structure
  cfg = cell2struct(val, var, 2);
end

% set the defaults
if ~isfield(cfg, 'location'),            cfg.location = 'interactive';     end
if ~isfield(cfg, 'locationcoordinates'), cfg.locationcoordinates = 'head'; end
if ~isfield(cfg, 'funparameter'),        cfg.funparameter = [];            end
if ~isfield(cfg, 'maskparameter'),       cfg.maskparameter = [];           end
if ~isfield(cfg, 'colormap'),            cfg.colormap = 'jet';             end
if ~isfield(cfg, 'downsample'),          cfg.downsample = 1;               end
if ~isfield(cfg, 'maskcolmin'),          cfg.maskcolmin = 'auto';          end
if ~isfield(cfg, 'maskcolmax'),          cfg.maskcolmax = 'auto';          end
if ~isfield(cfg, 'colmin'),              cfg.colmin = 'auto';              end
if ~isfield(cfg, 'colmax'),              cfg.colmax = 'auto';              end
if ~isfield(cfg, 'axis'),                cfg.axis   = 'on';                end
if ~isfield(cfg, 'crosshair'),           cfg.crosshair = 'on';             end
if ~isfield(cfg, 'colorbar'),            cfg.colorbar  = 'on';             end

% experimental options
if ~isfield(cfg, 'TTlookup');            cfg.TTlookup = 'no';              end
if ~isfield(cfg, 'TTqueryrange');        cfg.TTqueryrange = 3;             end

if isfield(cfg, 'flipdim')
  % FIXME the flipdim should be implemented cleanly, i.e. the volume
  % including the transformation should be flipped. For the moment I have
  % removed the flipdim code from this function.
  error('cfg.flipdim invalidates the location and voxel coordinates that are printed on screen\n');
end

if strcmp(cfg.TTlookup, 'yes')
  fprintf('reading Talairach-Tournoux atlas coordinates and labels\n');
  tlrc = TTatlas_init;
end

if isstr(interp)
  % read the anatomical MRI data from file
  filename = interp;
  fprintf('reading MRI from file\n');
  interp = read_fcdc_mri(filename);
end

% collect the volumes that will be plotted
ana = [];
fun = [];
msk = [];

% convert the coordinates along the axes (i.e. xgrid/ygrid/zgrid) into a homogenous transformation matrix
interp = grid2transform(interp);

% select the functional and the mask parameter
cfg.funparameter  = parameterselection(cfg.funparameter, interp);
cfg.maskparameter = parameterselection(cfg.maskparameter, interp);
% only a single parameter should be selected
try, cfg.funparameter  = cfg.funparameter{1};  end
try, cfg.maskparameter = cfg.maskparameter{1}; end

% downsample all volumes
tmpcfg = [];
tmpcfg.parameter  = {cfg.funparameter, cfg.maskparameter, 'anatomy'};
tmpcfg.downsample = cfg.downsample;
interp = volumedownsample(tmpcfg, interp);
interp.xgrid = 1:interp.dim(1);
interp.ygrid = 1:interp.dim(2);
interp.zgrid = 1:interp.dim(3);

if isfield(interp, 'anatomy')
  hasana = 1;
  mri8  = isa(interp.anatomy, 'uint8');
  mri16 = isa(interp.anatomy, 'uint16');
  % convert integers to single precision float if neccessary
  if mri8 || mri16
    fprintf('converting anatomy to double\n');
    ana = double(interp.anatomy);
  else
    ana = interp.anatomy;
  end
else
  hasana = 0;
  fprintf('no anatomical volume present\n');
  ana = ones(interp.dim);
end

if ~isempty(cfg.funparameter) && issubfield(interp, cfg.funparameter)
  hasfun = 1;
  fun = getsubfield(interp, cfg.funparameter);
else
  hasfun = 0;
  cfg.funparameter = [];
  fun = zeros(interp.dim);
  fprintf('no functional parameter\n');
end

% scale functional data
fprintf('scaling functional data...');
fmin = min(fun(:));
fmax = max(fun(:));
if ~ischar(cfg.colmin)
  fcolmin = cfg.colmin;
else
  if sign(fmin)==sign(fmax)
    fcolmin = fmin;
  else
    fcolmin = -max(abs([fmin,fmax]));
  end
end
if ~ischar(cfg.colmax)
  fcolmax = cfg.colmax;
else
  if sign(fmin)==sign(fmax)
    fcolmax = fmax;
  else
    fcolmax = max(abs([fmin,fmax]));
  end
end
fun = (fun-fcolmin)./(fcolmax-fcolmin);
if ~ischar(cfg.colmax)
  fun(find(fun>1)) = 1;
end
if ~ischar(cfg.colmin)
  fun(find(fun<0)) = 0;
end
fprintf('done\n');

if ~isempty(cfg.maskparameter) && issubfield(interp, cfg.maskparameter)
  hasmsk = 1;
  msk = getsubfield(interp, cfg.maskparameter);
  fprintf('scaling mask data...');
  mmin = min(msk(:));
  mmax = max(msk(:));
  if ~ischar(cfg.maskcolmin)
    mcolmin = cfg.maskcolmin;
  else
    if sign(mmin)==sign(mmax)
      mcolmin = mmin;
    else
      mcolmin = -max(abs([mmin,mmax]));
    end
  end
  if ~ischar(cfg.maskcolmax)
    mcolmax = cfg.maskcolmax;
  else
    if sign(mmin)==sign(mmax)
      mcolmax = mmax;
    else
      mcolmax = max(abs([mmin,mmax]));
    end
  end
  msk = (msk-mcolmin)./(mcolmax-mcolmin);
  if ~ischar(cfg.maskcolmax)
    msk(find(msk>1)) = 1;
  end
  if ~ischar(cfg.maskcolmin)
    msk(find(msk<0)) = 0;
  end
  fprintf('done\n');
else
  hasmsk = 0;
  cfg.maskparameter = [];
  if isempty(cfg.funparameter)
    msk = zeros(interp.dim);
  else
    msk = 0.5 * ones(interp.dim);
  end
  fprintf('no masking parameter\n');
end

% ensure that they are all 3D volumes
ana = reshape(ana, interp.dim);
fun = reshape(fun, interp.dim);
msk = reshape(msk, interp.dim);

% ensure that the functional data is real
if ~isreal(fun)
  fprintf('taking absolute value of complex data\n');
  fun = abs(fun);
end

anasc = [];
funsc = [];
msksc = [];

% automatically determine the interesting range of values that should be plotted
if isempty(anasc) && hasana
  anasc(1) = 0;
  anasc(2) = max(ana(:));
else
  anasc = [0 1];
end
if isempty(funsc) && hasfun
  funsc(1) = min(fun(:));
  funsc(2) = max(fun(:));
else
  funsc = [0 1];
end
if isempty(msksc) && hasmsk
  msksc(1) = 0;
  msksc(2) = 1;
else
  msksc = [0 1];
end

x = interp.xgrid;
y = interp.ygrid;
z = interp.zgrid;
dim = interp.dim;

if ~isstr(cfg.location)
  if strcmp(cfg.locationcoordinates, 'head')
    % convert the headcoordinates location into voxel coordinates
    sel = inv(interp.transform) * [cfg.location(:); 1];
    sel = round(sel(1:3));
  elseif strcmp(cfg.locationcoordinates, 'voxel')
    % the location is already in voxel coordinates
    sel = round(cfg.location(1:3));
  end
else
  sel = cfg.location;
end

% determine the initial intersection of the cursor
if isstr(sel) && strcmp(sel, 'min')
  if isempty(cfg.funparameter)
    error('no functional parameter specified');
  end
  [minval, minindx] = min(fun(:));
  [xi, yi, zi] = ind2sub(interp.dim, minindx);
elseif isstr(sel) && strcmp(sel, 'max')
  if isempty(cfg.funparameter)
    error('no functional parameter specified');
  end
  [maxval, maxindx] = max(fun(:));
  [xi, yi, zi] = ind2sub(interp.dim, maxindx);
elseif isstr(sel) && strcmp(sel, 'maxproject')
  if isempty(cfg.funparameter)
    error('no functional parameter specified');
  end
  % take [1 1 1] and project the maxima and contour of brain onto fun and ana.
  xi = 1;
  yi = 1;
  zi = 1; 
  
  load /home/coherence/jansch/matlab/fieldtrip/private/brainedges
 
  %create dummy anatomy with the outline of the brain in the first slices 
  ana        = zeros(size(ana));
  %postprocess edges
  sagittal(find(isnan(sagittal))) = 0;
  tmp                             = imfill(imdilate(sagittal==1, strel('diamond', 2)), [100 100]);
  sagittal                        = imerode(tmp, strel('diamond', 2));  
  ana(1,:,:)                      = double(~(tmp - sagittal));
  
  coronal(find(isnan(coronal))) = 0;
  tmp                           = imfill(imdilate(coronal==1, strel('diamond', 2)), [100 100]);
  coronal                       = imerode(tmp, strel('diamond', 2));  
  ana(:,1,:)                    = double(~(tmp - coronal));
  
  axial(find(isnan(axial))) = 0;
  tmp                       = imfill(imdilate(axial==1, strel('diamond', 2)), [100 100]);
  axial                     = imerode(tmp, strel('diamond', 2));  
  ana(:,:,1)                = double(~(tmp - axial));

  anasc = [0 1];
  
  %postprocess functional data and keep edges visible
  fun(1,:,:) = squeeze(max(fun, [], 1)).*double(sagittal);
  fun(:,1,:) = squeeze(max(fun, [], 2)).*double(coronal);
  fun(:,:,1) = squeeze(max(fun, [], 3)).*double(axial);
  
  msk(1,:,:) = squeeze(max(msk, [], 1)).*double(sagittal);
  msk(:,1,:) = squeeze(max(msk, [], 2)).*double(coronal);
  msk(:,:,1) = squeeze(max(msk, [], 3)).*double(axial);

elseif isstr(sel) && strcmp(sel, 'center')
  xi = round(length(x)/2);
  yi = round(length(y)/2);
  zi = round(length(z)/2);
elseif isstr(sel) && strcmp(sel, 'interactive')
  % start at the center
  xi = round(length(x)/2);
  yi = round(length(y)/2);
  zi = round(length(z)/2);
elseif ~isstr(sel)
  xi = nearest(x, sel(1));
  yi = nearest(y, sel(2));
  zi = nearest(z, sel(3));
end

% define the colormap for the functional data
funcolormap = colormap(cfg.colormap);

nas = [];
lpa = [];
rpa = [];

%if ~strcmp(cfg.location, 'interactive') && ~strcmp(cfg.location, 'maxproject'),
if ~strcmp(cfg.location, 'interactive'),
  % plot only once
  xi = round(xi); xi = max(xi, 1); xi = min(xi, dim(1));
  yi = round(yi); yi = max(yi, 1); yi = min(yi, dim(2));
  zi = round(zi); zi = max(zi, 1); zi = min(zi, dim(3));

  ijk = [xi yi zi 1]';
  xyz = interp.transform * ijk;
  val = fun(xi, yi, zi);
  fprintf('voxel %d, indices [%d %d %d], location [%.1f %.1f %.1f], value %f\n', sub2ind(interp.dim, xi, yi, zi), ijk(1:3), xyz(1:3), val);

  if strcmp(cfg.TTlookup, 'yes')
    lab = TTatlas_lookup(tlrc, mni2tal(xyz(1:3)), cfg.TTqueryrange);
    if isempty(lab)
      fprintf('Talairach-Tournoux labels: not found\n');
    else
      fprintf('Talairach-Tournoux labels: ')
      fprintf('%s', lab{1});
      for i=2:length(lab)
        fprintf(', %s', lab{i});
      end
      fprintf('\n');
    end
  end
  
  subplot(2,2,1); singleplot(ana, fun, msk, anasc, funsc, msksc, [xi yi zi], 2, funcolormap); xlabel('i'); ylabel('k'); axis(cfg.axis);
  if strcmp(cfg.crosshair, 'on'), crosshair([xi zi]); end
  subplot(2,2,2); singleplot(ana, fun, msk, anasc, funsc, msksc, [xi yi zi], 1, funcolormap); xlabel('j'); ylabel('k'); axis(cfg.axis);
  if strcmp(cfg.crosshair, 'on'), crosshair([yi zi]); end
  subplot(2,2,3); singleplot(ana, fun, msk, anasc, funsc, msksc, [xi yi zi], 3, funcolormap); xlabel('i'); ylabel('j'); axis(cfg.axis);
  if strcmp(cfg.crosshair, 'on'), crosshair([xi yi]); end
  if strcmp(cfg.colorbar,  'on'),
    vectorcolorbar = linspace(funsc(1),funsc(2),length(funcolormap));
    subplot(2,2,4);imagesc(vectorcolorbar,1,vectorcolorbar);colormap(funcolormap);
  end

  drawnow;
%elseif strcmp(cfg.location, 'maxproject'),
%  % make a 'glass-brain' projection
%  % this part is rough and experimental
%  xi = 1;
%  yi = 1;
%  zi = 1;
%  keyboard
%  load /home/coherence/jansch/matlab/fieldtrip/private/brainedges
%  ana(1,:,:) = sagittal;
%  ana(:,1,:) = coronal;
%  ana(:,:,1) = axial;  
%
%  % make a 'glass-brain' projection
%  % this part is rough and experimental
%  load /home/coherence/jansch/matlab/fieldtrip/private/brainedges
%  
%  %sagittal
%  [Xs,Ys]           = ndgrid(1:size(sagittal,1),1:size(sagittal,2));
%  Zs                = squeeze(max(fun, [], 1));
%  tmp               = sagittal;
%  tmp(find(tmp==0)) = nan;
%  tmp(find(~isnan(tmp)))    = 1;
%  subplot(2,2,2); hold on; contour(Xs, Ys, sagittal, 1, 'k', 'LineWidth', 2); axis equal  
%  subplot(2,2,2); h1 = surface(Xs, Ys, zeros(size(Zs)), Zs.*tmp, 'EdgeColor', 'none');
%  set(h1, 'AlphaData', squeeze(max(msk, [], 1)));
%  set(h1, 'FaceAlpha', 'interp');
%  axis off
% 
%  %coronal
%  [Xs,Ys]           = ndgrid(1:size(coronal,1),1:size(coronal,2));
%  Zs                = squeeze(max(fun, [], 2));
%  tmp               = coronal;
%  tmp(find(tmp==0)) = nan;
%  tmp(find(~isnan(tmp)))    = 1;
%  subplot(2,2,1); hold on; contour(Xs, Ys, coronal, 1, 'k', 'LineWidth', 2); axis equal  
%  subplot(2,2,1); h2 = surface(Xs, Ys, zeros(size(Zs)), Zs.*tmp, 'EdgeColor', 'none');
%  set(h2, 'AlphaData', msk);
%  set(h2, 'AlphaData', squeeze(max(msk, [], 2)));
%  set(h2, 'FaceAlpha', 'interp');
%  axis off  
%  
%  %axial
%  [Xs,Ys]           = ndgrid(1:size(axial,1),1:size(axial,2));
%  Zs                = squeeze(max(fun, [], 3));
%  tmp               = axial;
%  tmp(find(tmp==0)) = nan;
%  tmp(find(~isnan(tmp)))    = 1;
%  subplot(2,2,3); hold on; contour(Xs, Ys, axial, 1, 'k', 'LineWidth', 2); axis equal  
%  subplot(2,2,3); h3 = surface(Xs, Ys, zeros(size(Zs)), Zs.*tmp, 'EdgeColor', 'none');
%  set(h3, 'AlphaData', squeeze(max(msk, [], 3)));
%  set(h3, 'FaceAlpha', 'interp');
%  axis off
else
  % keep on plotting until the user presses quit
  while(1)

    xi = round(xi); xi = max(xi, 1); xi = min(xi, dim(1));
    yi = round(yi); yi = max(yi, 1); yi = min(yi, dim(2));
    zi = round(zi); zi = max(zi, 1); zi = min(zi, dim(3));

    fprintf('\n');
    fprintf('click with mouse button to reposition the cursor\n');
    fprintf('press q on keyboard to quit interactive mode\n');
    if ~isempty(nas), fprintf('nas = [%f %f %f]\n', nas); end
    if ~isempty(lpa), fprintf('lpa = [%f %f %f]\n', lpa); end
    if ~isempty(rpa), fprintf('rpa = [%f %f %f]\n', rpa); end

    ijk = [xi yi zi 1]';
    xyz = interp.transform * ijk;
    val = fun(xi, yi, zi);
    fprintf('voxel %d, indices [%d %d %d], location [%.1f %.1f %.1f], value %f\n', sub2ind(interp.dim, xi, yi, zi), ijk(1:3), xyz(1:3), val);

    if strcmp(cfg.TTlookup, 'yes')
      lab = TTatlas_lookup(tlrc, mni2tal(xyz(1:3)), cfg.TTqueryrange);
      if isempty(lab)
        fprintf('Talairach-Tournoux labels: not found\n');
      else
        fprintf('Talairach-Tournoux labels: ')
        fprintf('%s', lab{1});
        for i=2:length(lab)
          fprintf(', %s', lab{i});
        end
        fprintf('\n');
      end
    end

    subplot(2,2,1); cla; singleplot(ana, fun, msk, anasc, funsc, msksc, [xi yi zi], 2, funcolormap); xlabel('i'); ylabel('k'); crosshair([xi zi]);
    subplot(2,2,2); cla; singleplot(ana, fun, msk, anasc, funsc, msksc, [xi yi zi], 1, funcolormap); xlabel('j'); ylabel('k'); crosshair([yi zi]);
    subplot(2,2,3); cla; singleplot(ana, fun, msk, anasc, funsc, msksc, [xi yi zi], 3, funcolormap); xlabel('i'); ylabel('j'); crosshair([xi yi]);
    vectorcolorbar = linspace(funsc(1),funsc(2),length(funcolormap));
    subplot(2,2,4);imagesc(vectorcolorbar,1,vectorcolorbar);colormap(funcolormap);
    drawnow;

    try, [d1, d2, key] = ginput(1); catch, key='q'; end
    if key=='q'
      break;
    elseif key=='l'
      lpa = [xi yi zi];
    elseif key=='r'
      rpa = [xi yi zi];
    elseif key=='n'
      nas = [xi yi zi];
    else
      % update the view to a new position
      l1 = get(get(gca, 'xlabel'), 'string');
      l2 = get(get(gca, 'ylabel'), 'string');
      switch l1,
        case 'i'
          xi = d1;
        case 'j'
          yi = d1;
        case 'k'
          zi = d1;
      end
      switch l2,
        case 'i'
          xi = d2;
        case 'j'
          yi = d2;
        case 'k'
          zi = d2;
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINGLEPLOT makes an overlay of 3D anatomical, functional and probability
% volumes. The three volumes must be scaled between 0 and 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function singleplot(ana, fun, msk, anasc, funsc, msksc, indx, dimension ,funcolormap);

% select the indices of the intersection
xi = indx(1);
yi = indx(2);
zi = indx(3);

% select the slice to plot
if dimension==1
  yi = 1:size(ana,2);
  zi = 1:size(ana,3);
elseif dimension==2
  xi = 1:size(ana,1);
  zi = 1:size(ana,3);
elseif dimension==3
  xi = 1:size(ana,1);
  yi = 1:size(ana,2);
end

% cut out the slice of interest
ana = squeeze(ana(xi,yi,zi));
fun = squeeze(fun(xi,yi,zi));
msk = squeeze(msk(xi,yi,zi));

% fprintf('scaling and clipping anatomy\n');
ana(find(ana(:)<anasc(1))) = anasc(1);      % clip to minimal interesting values
ana(find(ana(:)>anasc(2))) = anasc(2);      % clip to maximal interesting values
ana = (ana-anasc(1))./(anasc(2)-anasc(1));  % scale interesting range to 0-1
% fprintf('scaling and clipping functional\n');
fun(find(fun(:)<funsc(1))) = funsc(1);      % clip to minimal values
fun(find(fun(:)>funsc(2))) = funsc(2);      % clip to maximal values
fun = (fun-funsc(1))./(funsc(2)-funsc(1));  % scale interesting range to 0-1
% fprintf('scaling and clipping mask\n');
msk(find(msk(:)<msksc(1))) = msksc(1);      % clip to minimal interesting values
msk(find(msk(:)>msksc(2))) = msksc(2);      % clip to maximal interesting values
msk = (msk-msksc(1))./(msksc(2)-msksc(1));  % scale interesting range to 0-1

% this is needed for displaying the matrix using the Matlab image() function
ana = ana';
fun = fun';
msk = msk';
dim = size(ana);

% convert anatomy into RGB values
ana = cat(3, ana, ana, ana);

% convert functional into RGB values
fun = floor((size(funcolormap,1)-1) * fun)+1;	% scale interesting range to 1-64
fun(find(isnan(fun(:)))) = 1;
r = zeros(dim);
g = zeros(dim);
b = zeros(dim);
r(:) = funcolormap(fun(:), 1);			% find the matching color (Rgb)
g(:) = funcolormap(fun(:), 2);			% find the matching color (rGb)
b(:) = funcolormap(fun(:), 3);			% find the matching color (rgB)
fun = cat(3, r, g, b);

ha = image(ana);
hold on
hf = image(fun);
set(hf, 'AlphaData', msk);          % apply the opacity mask to the functional data
axis equal
axis tight
axis xy


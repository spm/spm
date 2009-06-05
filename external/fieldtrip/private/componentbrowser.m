function [varargout] = componentbrowser(cfg, comp)

% COMPONENTBROWSER plots topography and activations of ICA components
%
% Use as
%   componentbrowser(cfg, comp)
% where comp is a FieldTrip structure obtained from COMPONENTANALYSIS.
%
% The configuration has the following parameters:
% cfg.layout   = layout from PREPARE_LAYOUT (required)
% cfg.comp  = a vector with the components to plot (ex. 1:10) (optional)
% cfg.trial = choose which trial to plot first (optional, only one trial)

% Copyright (C) 2009, Giovanni Piantoni
%
% $Log: componentbrowser.m,v $
% Revision 1.2  2009/06/03 14:00:26  roboos
% fixed cfg.lay, should be cfg.layout
%
% Revision 1.1  2009/06/02 15:48:58  giopia
% first implementation, plot topoplot, activations and simple interactive
%

fieldtripdefs

% check the data comes from componentanalysis
comp = checkdata(comp, 'datatype', 'comp');

% set the defaults:
if ~isfield(cfg, 'comp'),   cfg.comp  = 1:10;   end
if ~isfield(cfg, 'trial'),  cfg.trial   = 1;    end

% Read or create the layout that will be used for plotting:
lay = prepare_layout(cfg, comp);
cfg.layout = lay;

if numel(cfg.trial) > 1,
  warning('only one trial can be plotted at the time');
  cfg.trial = cfg.trial(1);
end

% fixed variables
cfg.shift   = 1.2;   % distance between topoplots
gridscale   = 67;   % default parameter from topoplot
[labels, chanidx] = intersect(comp.topolabel, lay.label); % in case channels are missing
[mask] = createmask(lay, gridscale);

% create figure
cfg.hndl = figure('uni','pix', 'name', 'componentbrowser', 'vis', 'off', 'numbertitle', 'off');
hold on

cnt = 0;
for k = cfg.comp
  cnt = cnt + 1;

  % write number of the component on the left
  plot_text(-2, -cnt*cfg.shift, ['n ' num2str(cfg.comp(cnt))])

  % plot only topography (no layout)
  plot_topo(comp.topo(chanidx, k), lay.pos(chanidx,1), lay.pos(chanidx,2), ...
    'hpos', -1, 'vpos', -cnt*cfg.shift, 'mask', mask, 'gridscale', gridscale, 'width', .5/.45, 'height', .5/.45) % .5/.45 is a magic number: necessary to adjust layout with mask
  % plot layout
  plot_lay(lay, 'hpos', -1, 'vpos', -cnt*cfg.shift, 'point', false, 'box', false, 'label', false, 'mask', true, 'verbose', false);
end

% plot activations
plotactivation([], cfg, comp) % first call of this function

% previous trial button
h_prev(1) = plot_box([5  7 -(cnt)*cfg.shift-1 -(cnt)*cfg.shift-1.5], 'facecolor', [.5 .5 .5]);
h_prev(2) = plot_text(6, -(cnt)*cfg.shift-1.4, '<<');
set(h_prev, 'ButtonDownFcn', {@plotactivation, comp}, 'tag', 'prev')

% next trial button
h_next(1) = plot_box([9 11 -(cnt)*cfg.shift-1 -(cnt)*cfg.shift-1.5], 'facecolor', [.5 .5 .5]);
h_next(2) = plot_text(10, -(cnt)*cfg.shift-1.4, '>>');
set(h_next, 'ButtonDownFcn', {@plotactivation, comp}, 'tag', 'next')

% final adjustments
set(cfg.hndl, 'vis', 'on')
axis equal
axis off
hold off

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = cfg.hndl;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTACTIVATION: subfunction which plots only the activation
% and it can be called more times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotactivation(varargin)
% plotactivation can be called in isolation or by buttondownfcn
% cfg is stored in 'user' of the main figure

comp    = varargin{3};

if isempty(varargin{1}) % when called in isolation
  cfg = varargin{2};
  set(cfg.hndl, 'user', cfg)
else
  cfg = get(get(get(varargin{1}, 'par'),'par'), 'user');

  % which button has been pressed
  if intersect(varargin{1}, findobj(cfg.hndl, 'tag', 'next'))
    cfg.trial = cfg.trial + 1;

    if cfg.trial > size(comp.trial,2)
      cfg.trial = size(comp.trial,2);
    end

  elseif intersect(varargin{1}, findobj(cfg.hndl, 'tag', 'prev'))
    cfg.trial = cfg.trial - 1;

    if cfg.trial < 1
      cfg.trial = 1;
    end

  end
end

delete(findobj(cfg.hndl,'tag', 'activations'));

hold on
cnt = 0;
for k = cfg.comp
  cnt = cnt + 1;

  % plot the activations
  h_act(cnt) = plot_vector(comp.trial{cfg.trial}(k,:), 'hpos', 6 , 'vpos', -cnt*cfg.shift, 'width', 12, 'height', 1, 'box', true);
end

% specify which trial is shown
h_act(cnt+1) = plot_text(8, -(cnt)*cfg.shift-1.4, ['trial n ' num2str(cfg.trial)]);

set(h_act, 'tag', 'activations')
set(cfg.hndl, 'user', cfg)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATEMASK: create anatomical mask, only one for all the topoplots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask] = createmask(lay, gridscale)
% calculate anatomical mask only once, based only on lay

% find limits for interpolation:
xmin = +inf;
xmax = -inf;
ymin = +inf;
ymax = -inf;
for i=1:length(lay.mask)
  xmin = min([xmin; lay.mask{i}(:,1)]);
  xmax = max([xmax; lay.mask{i}(:,1)]);
  ymin = min([ymin; lay.mask{i}(:,2)]);
  ymax = max([ymax; lay.mask{i}(:,2)]);
end

xi = linspace(xmin, xmax, gridscale);   % x-axis for interpolation (row vector)
yi = linspace(ymin, ymax, gridscale);   % y-axis for interpolation (row vector)
Xi =  ones(gridscale,1)*xi;
Yi = (ones(gridscale,1)*yi)';

% apply anatomical mask to the data, i.e. that determines that the interpolated data outside the circle is not displayed
mask = false(gridscale);
for i=1:length(lay.mask)
  lay.mask{i}(end+1,:) = lay.mask{i}(1,:); % force them to be closed
  mask(inside_contour([Xi(:) Yi(:)], lay.mask{i})) = true;
end



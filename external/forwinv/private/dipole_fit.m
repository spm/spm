function [dipout] = dipole_fit(dip, sens, vol, dat, varargin)

% DIPOLE_FIT performs an equivalent current dipole fit with a single
% or a small number of dipoles to explain an EEG or MEG scalp topography.
%
% Use as
%   [dipout] = dipole_fit(dip, sens, vol, dat, ...)
%
% Additional input arguments should be specified as key-value pairs and can include
%   'constr'      = Structure with constraints
%   'display'     = Level of display [ off | iter | notify | final ]
%   'optimfun'    = Function to use [fminsearch | fminunc ]
%   'maxiter'     = Maximum number of function evaluations allowed [ positive integer ]
%   'metric'      = Error measure to be minimised [ rv | var | abs ]
%   'checkinside' = Boolean flag to check whether dipole is inside source compartment [ 0 | 1 ]
%   'weight'      = weight matrix for maximum likelihood estimation, e.g. inverse noise covariance
%
% The following optional input arguments relate to the computation of the leadfields
%   'reducerank'      = 'no' or number
%   'normalize'       = 'no', 'yes' or 'column'
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%
% The maximum likelihood estimation implements
%   Lutkenhoner B. "Dipole source localization by means of maximum
%   likelihood estimation I. Theory and simulations" Electroencephalogr Clin
%   Neurophysiol. 1998 Apr;106(4):314-21.

% Copyright (C) 2003-2008, Robert Oostenveld
%
% $Log: dipole_fit.m,v $
% Revision 1.11  2008/04/30 13:48:57  roboos
% finished implementation of MLE fitting, renamed cov into weight (inverse should be done outside this function)
%
% Revision 1.10  2008/04/29 15:41:53  roboos
% added skelleton for MAP estimation, but the covariance is not passed yet to lower level function
%
% Revision 1.9  2008/01/09 12:43:25  roboos
% added optional input arguments for reducerank, normalize and normalizeparam, these are directly passed on to compute_leadfield
%
% Revision 1.8  2007/12/18 13:01:55  roboos
% added support for fitting fixed orientation dipoles (constr.fixedori=1)
%
% Revision 1.7  2007/07/23 10:16:00  roboos
% implemented options for changing the error metric and for checking taht the dipole is inside the brain
%
% Revision 1.6  2007/03/21 11:13:29  roboos
% added MaxFunEvals to the optimization options
%
% Revision 1.5  2007/01/02 10:51:23  roboos
% renamed the output argument to dipout, no functional changes
%
% Revision 1.4  2006/09/07 12:39:51  roboos
% included optional input arguments to documentation
%
% Revision 1.3  2006/05/11 07:16:44  roboos
% Added additional input argument 'optimfun', can be fminunc or
% fminsearch. Changed the default maxiter to 500 for fminsearch.
% Changed the assignment of sens.pnt as structure (prevents warning
% in matlab7).
%
% Revision 1.2  2006/05/10 10:48:21  roboos
% Changed the order of the input arguments for consistency with
% other inverse functions, and switched to specifying optional inputs
% using key-value pairs. The old-style input is supported for backward
% compatibility with EEGLAB and dipfit1. Added optional input arguments
% for maxiter (default=100, just as before) and display (default='iter',
% just as before).
%
% Revision 1.1  2005/11/01 09:34:27  roboos
% replaced EEG/MEG specific dipole fit routines with a single general function
%
% Revision 1.7  2004/02/05 11:55:30  roberto
% fixed clearly apparent bug in error function (would not run)
%
% Revision 1.6  2003/09/12 08:45:41  roberto
% added error if maximum number of iterations is exceeded
%
% Revision 1.5  2003/09/02 13:01:18  roberto
% implemented constrained dipole fitting
%
% Revision 1.4  2003/06/16 10:03:31  roberto
% added check in error function for gui interrupt request (for eeglab)
%
% Revision 1.3  2003/03/13 13:41:00  roberto
% fixed bug in e/meg_error_func in assignment of Nchan
%
% Revision 1.2  2003/03/11 14:45:36  roberto
% updated help and copyrights
%

% It is neccessary to provide backward compatibility support for the old function call
% in case people want to use it in conjunction with EEGLAB and the dipfit1 plugin.
% old style: function [dipout] = dipole_fit(dip, dat, sens, vol, constr), where constr is optional
% new style: function [dipout] = dipole_fit(dip, sens, vol, dat, varargin), where varargin is in key-value pairs
if nargin==4 && ~isstruct(sens) && isstruct(dat)
  % looks like old style, the order of the input arguments has to be changed
  warning('converting from old style input\n');
  olddat   = sens;
  oldsens  = vol;
  oldvol   = dat;
  dat      = olddat;
  sens     = oldsens;
  vol      = oldvol;
elseif nargin==5  && ~isstruct(sens) && isstruct(dat)
  % looks like old style, the order of the input arguments has to be changed
  % furthermore the additional constraint has to be fixed
  warning('converting from old style input\n');
  olddat   = sens;
  oldsens  = vol;
  oldvol   = dat;
  dat      = olddat;
  sens     = oldsens;
  vol      = oldvol;
  varargin = {'constr', varargin{1}};  % convert into a  key-value pair
else
  % looks like new style, i.e. with optional key-value arguments
  % this is dealt with below
end

constr         = keyval('constr',         varargin); % default is not to have constraints
metric         = keyval('metric',         varargin); if isempty(metric),      metric = 'rv';     end
checkinside    = keyval('checkinside',    varargin); if isempty(checkinside), checkinside = 0;   end
display        = keyval('display',        varargin); if isempty(display),     display = 'iter';  end
optimfun       = keyval('optimfun',       varargin);
maxiter        = keyval('maxiter',        varargin);
reducerank     = keyval('reducerank',     varargin); % for leadfield computation
normalize      = keyval('normalize' ,     varargin); % for leadfield computation
normalizeparam = keyval('normalizeparam', varargin); % for leadfield computation
weight            = keyval('weight',            varargin); % for maximum likelihood estimation

% determine whether the Matlab Optimization toolbox can be used
if isempty(optimfun)
  if exist('fminunc')
    optimfun = 'fminunc';
  else
    optimfun = 'fminsearch';
  end
end

% set a default for the maximum number of iterations, depends on the optimization function
if isempty(maxiter)
  if strcmp(optimfun, 'fminunc')
    maxiter = 100;
  else
    maxiter = 500;
  end
end

% determine whether it is EEG or MEG
if ~isstruct(sens) & size(sens,2)==3
  % definition of electrode positions only, restructure it
  sens = struct('pnt', sens);
  iseeg = 1;
  ismeg = 0;
elseif isfield(sens, 'type') & strcmp(sens.type, 'electrodes')
  iseeg = 1;
  ismeg = 0;
elseif isfield(sens, 'type') & strcmp(sens.type, 'gradiometers')
  iseeg = 0;
  ismeg = 1;
elseif isfield(sens, 'pnt') & ~isfield(sens, 'ori')
  iseeg = 1;
  ismeg = 0;
elseif isfield(sens, 'pnt') &  isfield(sens, 'ori')
  iseeg = 0;
  ismeg = 1;
else
  iseeg = 0;
  ismeg = 0;
end

if ismeg & iseeg
  % this is something that I might implement in the future
  error('simultaneous EEG and MEG not supported');
elseif iseeg
  % ensure that the potential data is average referenced, just like the model potential
  dat = avgref(dat);
end

% reformat the position parameters in case of multiple dipoles
param = dip.pos(:)';
numdip = numel(param)/3;

% add the orientation to the nonlinear parameters
if isfield(constr, 'fixedori') && constr.fixedori
  for i=1:numdip
    % add the orientation to the list of parameters
    [th, phi, r] = cart2sph(dip.mom(1,i), dip.mom(2,i), dip.mom(3,i));
    param = [param th phi];
  end
end

% reduce the number of parameters to be fitted according to the constraints
if isfield(constr, 'mirror')
  param = param(constr.reduce);
end

if strcmp(optimfun, 'fminunc')
  % do non-linear optimization of the dipole parameters using Optimization toolbox fminunc()
  options = optimset('TolFun',1e-9,...
    'TypicalX',ones(size(param)),...
    'LargeScale','off',...
    'HessUpdate','bfgs',...
    'MaxIter',maxiter,...
    'MaxFunEvals',2*maxiter*length(param),...
    'Display',display);
  [param, fval, exitflag, output] = fminunc(@dipfit_error, param, options, dat, sens, vol, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight);
else
  % do non-linear optimization of the dipole parameters using standard Matlab fminsearch()
  options = optimset('MaxIter',maxiter,...
    'MaxFunEvals',2*maxiter*length(param),...
    'Display',display);
  [param, fval, exitflag, output] = fminsearch(@dipfit_error, param, options, dat, sens, vol, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight);
end

if exitflag==0
  error('Maximum number of iterations exceeded before reaching the minimum, I suggest that you try with another initial guess.')
end

% do linear optimization of dipole moment parameters
[err, mom] = dipfit_error(param, dat, sens, vol, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight);

% expand the number of parameters according to the constraints
if isfield(constr, 'mirror')
  param = constr.mirror .* param(constr.expand);
end

% get the dipole position and orientation
if isfield(constr, 'fixedori') && constr.fixedori
  numdip = numel(param)/5;
  ori = zeros(3,numdip);
  for i=1:numdip
    th  = param(end-(2*i)+1);
    phi = param(end-(2*i)+2);
    [ori(1,i), ori(2,i), ori(3,i)] = sph2cart(th, phi, 1);
  end
  pos = reshape(param(1:(numdip*3)), 3, numdip)';
else
  numdip = numel(param)/3;
  pos = reshape(param, 3, numdip)';
end

% return the optimal dipole parameters
dipout.pos = pos;
if isfield(constr, 'fixedori') && constr.fixedori
  dipout.mom  = ori;  % dipole orientation as vector
  dipout.ampl = mom;  % dipole strength
else
  dipout.mom  = mom;  % dipole moment, which represents both the orientation and strength as vector
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIPFIT_ERROR computes the error between measured and model data
% and can be used for non-linear fitting of dipole position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, mom] = dipfit_error(param, dat, sens, vol, constr, metric, checkinside, reducerank, normalize, normalizeparam, weight)

% flush pending graphics events, ensure that fitting is interruptible
drawnow;
if ~isempty(get(0, 'currentfigure')) & strcmp(get(gcf, 'tag'), 'stop')
  % interrupt the fitting
  close;
  error('USER ABORT');
end;

% expand the number of parameters according to the constraints
if isfield(constr, 'mirror')
  param = constr.mirror .* param(constr.expand);
end

% get the dipole positions and optionally also the orientation
if isfield(constr, 'fixedori') && constr.fixedori
  numdip = numel(param)/5;
  ori = zeros(3,numdip);
  for i=1:numdip
    th  = param(end-(2*i)+1);
    phi = param(end-(2*i)+2);
    [ori(1,i), ori(2,i), ori(3,i)] = sph2cart(th, phi, 1);
  end
  pos = reshape(param(1:(numdip*3)), 3, numdip)';
else
  numdip = numel(param)/3;
  pos = reshape(param, 3, numdip)';
end

% check whether the dipole is inside the source compartment
if checkinside
  [inside, outside] = find_inside_vol(pos, vol);
  if ~isempty(outside)
    error('Dipole is outside the source compartment');
  end
end

% construct the leadfield matrix for all dipoles
lf = compute_leadfield(pos, sens, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
if isfield(constr, 'fixedori') && constr.fixedori
  lf = lf * ori;
end

% compute the optimal dipole moment and the model error
if ~isempty(weight)
  % maximum likelihood estimation using the weigth matrix
  mom = pinv(lf'*weight*lf)*lf'*weight*dat;  % Lutkenhoner equation 5
  dif = dat - lf*mom;
  % compute the generalized goodness-of-fit measure
  switch metric
    case 'rv' % relative residual variance
      num   = dif' * weight * dif;
      denom = dat' * weight * dat;
      err   = sum(num(:)) ./ sum(denom(:)); % Lutkenhonner equation 7, except for the gof=1-rv
    case 'var' % residual variance
      num   = dif' * weight * dif';
      err   = sum(num(:));
    otherwise
      error('Unsupported error metric for maximum likelihood dipole fitting');
  end
else
  % ordinary least squares, this is the same as MLE with weight=eye(nchans,nchans)
  mom = pinv(lf)*dat;
  dif = dat - lf*mom;
  % compute the ordinary goodness-of-fit measures
  switch metric
    case 'rv' % relative residual variance
      err = sum(dif(:).^2) / sum(dat(:).^2);
    case 'var' % residual variance
      err = sum(dif(:).^2);
    case 'abs' % absolute difference
      err = sum(abs(dif));
    otherwise
      error('Unsupported error metric for dipole fitting');
  end
end


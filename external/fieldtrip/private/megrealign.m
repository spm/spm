function [interp] = megrealign(cfg, data);

% MEGREALIGN interpolates MEG data towards standard gradiometer locations 
% by projecting the individual timelocked data towards a coarse source
% reconstructed representation and computing the magnetic field on
% the standard gradiometer locations.
%
% Use as
%   [interp] = megrealign(cfg, data)
%   Required configuration options:
%   cfg.template, cfg.inwardshift
%
% The new gradiometer definition is obtained from a template dataset,
% or can be constructed by averaging the gradiometer positions over
% multiple datasets.
%   cfg.template       = single dataset that serves as template
%   cfg.template(1..N) = datasets that are averaged into the standard
%
% The realignment is done by computing a minumum current estimate using a
% large number of dipoles that are placed in the upper layer of the brain
% surface, followed by a forward computation towards the template
% gradiometer array. This requires the specification of a volume conduction
% model of the head and of a source model.
%
% A head model must be specified with
%   cfg.hdmfile     = string, file containing the volume conduction model
% or alternatively manually using 
%   cfg.vol.r       = radius of sphere
%   cfg.vol.o       = [x, y, z] position of origin
%
% A source model (i.e. a superficial layer with distributed sources) can be
% constructed from a headshape file, or from the volume conduction model
%   cfg.headshape   = filename for headshape, can be empty (default = [])
%   cfg.spheremesh  = number of dipoles in the source layer (default = 642)
%   cfg.inwardshift = depth of the source layer relative to the headshape 
%                     surface or volume conduction model (no default 
%                     supplied, see below)
%
% If you specify a headshape file and it contains a skin surface, the
% inward shift should be 2.5.
%
% For a single-sphere or a local-spheres headmodel based on the skin
% surface, an inward shift of 2.5 is reasonable. 
% 
% For a single-sphere or a local-spheres headmodel based on the brain
% surface, you should probably use an inward shift of about 1.
% 
% For a realistic single-shell headmodel based on the brain surface, you
% should probably use an inward shift of about 1. 
% 
% Other options are
%   cfg.pruneratio  = for singular values, default is 1e-3
%   cfg.verify      = 'yes' or 'no', show the percentage difference (default = 'yes')
%   cfg.feedback    = 'yes' or 'no' (default = 'no')
%   cfg.channel     =  Nx1 cell-array with selection of channels (default = 'MEG'),
%                      see CHANNELSELECTION for details
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
% 
% This implements the method described by T.R. Knosche, Transformation
% of whole-head MEG recordings between different sensor positions.
% Biomed Tech (Berl). 2002 Mar;47(3):59-62.
% 
% See also MEGINTERPOLATE, PREPARE_LOCALSPHERES, PREPARE_SINGLESHELL

% This function depends on PREPARE_DIPOLE_GRID
%
% This function depends on PREPARE_VOL_SENS which has the following options:
% cfg.channel
% cfg.elec
% cfg.elecfile
% cfg.grad
% cfg.gradfile
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented

% Copyright (C) 2004-2007, Robert Oostenveld
%
% $Log: megrealign.m,v $
% Revision 1.56  2008/12/02 12:22:00  roboos
% allow data to be realigned AND simultaneously interpolated to another MEG sensor type, e.g. from ctf151 to ctf275
%
% Revision 1.55  2008/11/25 14:56:52  estmee
% Documentation update
%
% Revision 1.54  2008/11/21 12:48:17  sashae
% added call to checkconfig at start and end of function
%
% Revision 1.53  2008/10/02 15:32:21  sashae
% replaced call to createsubcfg with checkconfig
%
% Revision 1.52  2008/09/30 16:45:55  sashae
% checkconfig: checks if the input cfg is valid for this function
%
% Revision 1.51  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.50  2008/07/15 19:56:44  roboos
% moved cfg details for dipole grid to subcfg (cfg.grid)subcfg (cfg.grid.xxx)

fieldtripdefs

cfg = checkconfig(cfg);

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'ismeg', 'yes');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'renamed',     {'plot3d',      'feedback'});
cfg = checkconfig(cfg, 'renamedval',  {'headshape',   'headmodel', []});
cfg = checkconfig(cfg, 'required',    {'inwardshift', 'template'});

% set the default configuration 
if ~isfield(cfg, 'headshape'),     cfg.headshape = [];            end
if ~isfield(cfg, 'pruneratio'),    cfg.pruneratio = 1e-3;         end
if ~isfield(cfg, 'spheremesh'),    cfg.spheremesh = 642;          end
if ~isfield(cfg, 'verify'),        cfg.verify = 'yes';            end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'yes';          end
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';            end
if ~isfield(cfg, 'channel'),       cfg.channel = 'MEG';           end
if ~isfield(cfg, 'topoparam'),     cfg.topoparam = 'rms';         end

% put the low-level options pertaining to the dipole grid in their own field
cfg = checkconfig(cfg, 'createsubcfg',  {'grid'});

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
  % update the trial definition (trl)
  if isfield(data, 'cfg') % try to locate the trl in the nested configuration
    trl = findcfg(data.cfg, 'trl');
  else
    trl = [];
  end
  if isempty(trl)
    % a trial definition is expected in each continuous data set
    warning('could not locate the trial definition ''trl'' in the data structure');
  else
    cfg.trlold=trl;
    cfg.trl=trl(cfg.trials,:);
  end
end

Ntrials = length(data.trial);

% retain only the MEG channels in the data and temporarily store
% the rest, these will be added back to the transformed data later.
cfg.channel = channelselection(cfg.channel, data.label);
dataindx = match_str(data.label, cfg.channel);
restindx = setdiff(1:length(data.label),dataindx);
if ~isempty(restindx)
  fprintf('removing %d non-MEG channels from the data\n', length(restindx));
  rest.label = data.label(restindx);    % first remember the rest
  data.label = data.label(dataindx);    % then reduce the data
  for i=1:Ntrials
    rest.trial{i} = data.trial{i}(restindx,:);	% first remember the rest
    data.trial{i} = data.trial{i}(dataindx,:);	% then reduce the data
  end
else
  rest.label = {};
  rest.trial = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the average template gradiometer array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ntemplate = length(cfg.template);
for i=1:Ntemplate
  fprintf('reading template helmet position from %s\n', cfg.template{i});
  template(i) = read_sens(cfg.template{i});
end

% to construct the average location of the MEG sensors, 4 channels are needed that should  be sufficiently far apart
switch sensortype(template(1))
case {'ctf151' 'ctf275'}
  labC = 'MZC01';
  labF = 'MZF03';
  labL = 'MLC21';
  labR = 'MRC21';
case {'ctf151_planar' 'ctf275_planar'}
  labC = 'MZC01_dH';
  labF = 'MZF03_dH';
  labL = 'MLC21_dH';
  labR = 'MRC21_dH';
case {'bti148' 'bti248'}
  labC = 'A14';
  labF = 'A2';
  labL = 'A15';
case {'bti148_planar' 'bti248_planar'}
  labC = 'A14';
  labF = 'A2';
  labL = 'A15';
  labR = 'A29';
otherwise
  % this could in principle be added to the cfg, but better is to have a more exhaustive list here
  error('unsupported MEG system for realigning, please ask on the mailing list');
end

templ.meanC = [0 0 0];
templ.meanF = [0 0 0];
templ.meanL = [0 0 0];
templ.meanR = [0 0 0];
for i=1:Ntemplate
  % determine the 4 ref sensors for this individual template helmet 
  indxC = strmatch(labC, template(i).label, 'exact'); 
  indxF = strmatch(labF, template(i).label, 'exact'); 
  indxL = strmatch(labL, template(i).label, 'exact'); 
  indxR = strmatch(labR, template(i).label, 'exact'); 
  if isempty(indxC) || isempty(indxF) || isempty(indxL) || isempty(indxR)
    error('not all 4 sensors were found that are needed to rotate/translate');
  end
  % add them to the sum, to compute mean location of each ref sensor
  templ.meanC = templ.meanC + template(i).pnt(indxC,:);
  templ.meanF = templ.meanF + template(i).pnt(indxF,:);
  templ.meanL = templ.meanL + template(i).pnt(indxL,:);
  templ.meanR = templ.meanR + template(i).pnt(indxR,:);
end
templ.meanC = templ.meanC / Ntemplate;
templ.meanF = templ.meanF / Ntemplate;
templ.meanL = templ.meanL / Ntemplate;
templ.meanR = templ.meanR / Ntemplate;

% construct two direction vectors that define the helmet orientation
templ.dirCF = (templ.meanF - templ.meanC);
templ.dirRL = (templ.meanL - templ.meanR);
% construct three orthonormal direction vectors
templ.dirX = normalize(templ.dirCF);
templ.dirY = normalize(templ.dirRL - dot(templ.dirRL, templ.dirX) * templ.dirX);
templ.dirZ = cross(templ.dirX, templ.dirY);
templ.tra = fixedbody(templ.meanC, templ.dirX, templ.dirY, templ.dirZ);

% determine the 4 ref sensors for the helmet that belongs to this dataset
indxC = strmatch(labC, data.grad.label, 'exact'); 
indxF = strmatch(labF, data.grad.label, 'exact'); 
indxL = strmatch(labL, data.grad.label, 'exact'); 
indxR = strmatch(labR, data.grad.label, 'exact'); 
if isempty(indxC) | isempty(indxF) | isempty(indxL) | isempty(indxR)
  error('not all 4 sensors were found that are needed to rotate/translate');
end

% construct two direction vectors that define the helmet orientation
orig.dirCF = template(1).pnt(indxF,:) - template(1).pnt(indxC,:);
orig.dirRL = template(1).pnt(indxL,:) - template(1).pnt(indxR,:);
% construct three orthonormal direction vectors
orig.dirX = normalize(orig.dirCF);
orig.dirY = normalize(orig.dirRL - dot(orig.dirRL, orig.dirX) * orig.dirX);
orig.dirZ = cross(orig.dirX, orig.dirY);
orig.tra = fixedbody(template(1).pnt(indxC,:), orig.dirX, orig.dirY, orig.dirZ);

% compute the homogenous transformation matrix and transform the positions
tra = inv(templ.tra) * orig.tra;
pnt = template(1).pnt;
pnt(:,4) = 1;
pnt = (tra * pnt')';
pnt = pnt(:,1:3);

% remove the translation from the transformation matrix and rotate the orientations
tra(:,4) = [0 0 0 1]';
ori = template(1).ori;
ori(:,4) = 1;
ori = (tra * ori')';
ori = ori(:,1:3);

tmp_label = template(1).label;
tmp_tra   = template(1).tra;
tmp_unit  = template(1).unit;

% construct the final template gradiometer definition
template = [];
template.grad.pnt = pnt;
template.grad.ori = ori;
% keep the same labels and the linear weights of the coils
template.grad.label = tmp_label;
template.grad.tra   = tmp_tra;
template.grad.unit  = tmp_unit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE_VOL_SENS will match the data labels, the gradiometer labels and the 
% volume model labels (in case of a multisphere model) and result in a gradiometer 
% definition that only contains the gradiometers that are present in the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpcfg = [];
if isfield(cfg, 'hdmfile')
  tmpcfg.hdmfile = cfg.hdmfile;
elseif isfield(cfg, 'vol') 
  tmpcfg.vol = cfg.vol;
end
tmpcfg.grad    = data.grad;
tmpcfg.channel = data.label; % this might be a subset of the MEG channels
[volold, data.grad] = prepare_headmodel(tmpcfg);

% note that it is neccessary to keep the two volume conduction models
% seperate, since the single-shell Nolte model contains gradiometer specific
% precomputed parameters. Note that this is not guaranteed to result in a
% good projection for local sphere models. 
tmpcfg.grad    = template.grad;
tmpcfg.channel = 'MEG'; % include all MEG channels
[volnew, template.grad] = prepare_headmodel(tmpcfg);

if strcmp(senstype(data.grad), senstype(template.grad))
  fprintf('mean distance towards template gradiometers is %.2f %s\n', mean(sum((data.grad.pnt-template.grad.pnt).^2, 2).^0.5), template.grad.unit);
else
  % the projection is from one MEG system to another MEG system, which makes a comparison of the data difficult
  cfg.feedback = 'no';
  cfg.verify = 'no';
end

% create the dipole grid on which the data will be projected
grid = prepare_dipole_grid(cfg, volold, data.grad);
pos = grid.pos;

% sometimes some of the dipole positions are nan, due to problems with the headsurface triangulation
% remove them to prevent problems with the forward computation
sel = find(any(isnan(pos(:,1)),2));
pos(sel,:) = [];

% compute the forward model for the old and new gradiometer positions
fprintf('computing forward model for %d dipoles\n', size(pos,1));
lfold = compute_leadfield(pos, data.grad,     volold);
lfnew = compute_leadfield(pos, template.grad, volnew);

% compute this inverse only once, although it is used twice
tmp = prunedinv(lfold, cfg.pruneratio);
% compute the three interpolation matrices
fprintf('computing interpolation matrix #1\n');
realign = lfnew * tmp;
if strcmp(cfg.verify, 'yes')
  fprintf('computing interpolation matrix #2\n');
  noalign = lfold * tmp;
  fprintf('computing interpolation matrix #3\n');
  bkalign = lfold * prunedinv(lfnew, cfg.pruneratio) * realign;
end

% interpolate the data towards the template gradiometers
for i=1:Ntrials
  fprintf('realigning trial %d\n', i);
  data.realign{i} = realign * data.trial{i};
  if strcmp(cfg.verify, 'yes')
    % also compute the residual variance when interpolating
    rvrealign = rv(data.trial{i}, data.realign{i});
    fprintf('original -> template             RV %.2f %%\n', 100 * mean(rvrealign));
    datnoalign = noalign * data.trial{i};
    datbkalign = bkalign * data.trial{i};
    rvnoalign = rv(data.trial{i}, datnoalign);
    rvbkalign = rv(data.trial{i}, datbkalign);
    fprintf('original             -> original RV %.2f %%\n', 100 * mean(rvnoalign));
    fprintf('original -> template -> original RV %.2f %%\n', 100 * mean(rvbkalign)); 
  end
end

% plot the topography before and after the realignment
if strcmp(cfg.feedback, 'yes')
    
  warning('showing MEG topography (RMS value over time) in the first trial only');
  Nchan = length(data.grad.label);
  pnt1 = data.grad.pnt(1:Nchan,:);
  pnt2 = template.grad.pnt(1:Nchan,:);
  prj1 = elproj(pnt1); tri1 = delaunay(prj1(:,1), prj1(:,2));
  prj2 = elproj(pnt2); tri2 = delaunay(prj2(:,1), prj2(:,2));

  switch cfg.topoparam
    case 'rms'
      p1 = sqrt(mean(data.trial{1}.^2, 2));
      p2 = sqrt(mean(data.realign{1}.^2, 2));
    case 'svd'
      [u, s, v] = svd(data.trial{1}); p1 = u(:,1);
      [u, s, v] = svd(data.realign{1}); p2 = u(:,1);
    otherwise
      error('unsupported cfg.topoparam');
  end

  X = [pnt1(:,1) pnt2(:,1)]';
  Y = [pnt1(:,2) pnt2(:,2)]';
  Z = [pnt1(:,3) pnt2(:,3)]';

  % show figure with old an new helmets, volume model and dipole grid
  figure
  tmpcfg = [];
  tmpcfg.vol = volold;
  tmpcfg.grad = data.grad;
  tmpcfg.grid = grid;
  tmpcfg.plotsensors = 'no';  % these are plotted seperately below
  headmodelplot(tmpcfg);
  hold on
  plot3(pnt1(:,1), pnt1(:,2), pnt1(:,3), 'r.') % original positions
  plot3(pnt2(:,1), pnt2(:,2), pnt2(:,3), 'g.') % template positions
  line(X,Y,Z, 'color', 'black');
  view(-90, 90);

  % show figure with data on old helmet location
  figure
  hold on
  plot3(pnt1(:,1), pnt1(:,2), pnt1(:,3), 'r.') % original positions
  plot3(pnt2(:,1), pnt2(:,2), pnt2(:,3), 'g.') % template positions
  line(X,Y,Z, 'color', 'black');
  axis equal; axis vis3d
  triplot(pnt1, tri1, p1);
  title('RMS, before realignment')
  view(-90, 90)
 
  % show figure with data on new helmet location
  figure
  hold on
  plot3(pnt1(:,1), pnt1(:,2), pnt1(:,3), 'r.') % original positions
  plot3(pnt2(:,1), pnt2(:,2), pnt2(:,3), 'g.') % template positions
  line(X,Y,Z, 'color', 'black');
  axis equal; axis vis3d
  triplot(pnt2, tri2, p2);
  title('RMS, after realignment')
  view(-90, 90)
end

% store the realigned data in a new structure
interp.label   = template.grad.label;
interp.grad    = template.grad;	  % replace with the template gradiometer array
interp.trial   = data.realign;    % remember the processed data
interp.fsample = data.fsample;
interp.time    = data.time;

% add the rest channels back to the data, these were not interpolated
if ~isempty(rest.label)
  fprintf('adding %d non-MEG channels back to the data (', length(rest.label));
  fprintf('%s, ', rest.label{1:end-1});
  fprintf('%s)\n', rest.label{end});
  for trial=1:length(rest.trial)
    interp.trial{trial} = [interp.trial{trial}; rest.trial{trial}];
  end
  interp.label = [interp.label; rest.label];
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% store the configuration of this function call, including that of the previous function call
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: megrealign.m,v 1.56 2008/12/02 12:22:00 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
interp.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction that computes the inverse using a pruned SVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lfi] = prunedinv(lf, r)
[u, s, v] = svd(lf);
p = find(s<(s(1,1)*r) & s~=0);
fprintf('pruning %d from %d, i.e. removing the %d smallest spatial components\n', length(p), min(size(s)), length(p));
s(p) = 0;
s(find(s~=0)) = 1./s(find(s~=0));
lfi = v * s' * u';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction that computes the homogenous translation matrix
% corresponding to a fixed body rotation and translation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = fixedbody(center, dirx, diry, dirz);
rot = eye(4);
rot(1:3,1:3) = inv(eye(3) / [dirx; diry; dirz]);
tra = eye(4);
tra(1:4,4)   = [-center 1]';
h = rot * tra;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction that scales a vector to unit length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v] = normalize(v);
v = v / sqrt(v * v');

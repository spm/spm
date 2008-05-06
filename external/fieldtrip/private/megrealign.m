function [interp] = megrealign(cfg, data);

% MEGREALIGN interpolates MEG data towards standard gradiometer locations 
% by projecting the individual timelocked data towards a coarse source
% reconstructed representation and computing the magnetic field on
% the standard gradiometer locations.
%
% Use as
%   [interp] = megrealign(cfg, data)
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
% Revision 1.49  2008/05/06 15:43:46  sashae
% change in trial selection, cfg.trials can be logical
%
% Revision 1.48  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.47  2008/03/05 10:46:36  roboos
% moved electrode reading functionality from read_fcdc_elec to read_sens, switched to the use of the new function
%
% Revision 1.46  2008/01/31 17:20:17  sashae
% added option for trial selection
%
% Revision 1.45  2007/10/15 10:53:21  roboos
% fixed bug in realignment when single-shell (a.k.a. "Nolte") headmodel
% was used. In the previous version, the data would not be realigned to
% the new template helmet, but would only be a little bit spatial noise
% filtered (i.e. slightly different).
%
% Revision 1.44  2007/05/30 13:23:55  roboos
% added cfg.channel default to ensure that only MEG channels (and not reference channels) are used for realignment. All other channels are not realigned.
%
% Revision 1.43  2007/05/29 16:10:24  ingnie
% use read_fcdc_elec iso read_fcdc_header to read gradiometer positions (with roboos)
%
% Revision 1.42  2007/05/29 12:51:31  roboos
% added new options for checkdata
%
% Revision 1.41  2007/05/16 12:22:04  roboos
% fixed some small bugs
%
% Revision 1.40  2007/05/16 11:46:20  roboos
% switched to using prepare_dipole_grid subfunction for source layer generation
% changed the plotting, an  additional figure is made with the geometry of the volume conductor, teh sources and the old+new sensor locations
% changed default cfg.feedback to yes (i.e. default is to show the figures)
%
% Revision 1.39  2007/05/02 15:59:13  roboos
% be more strict on the input and output data: It is now the task of
% the private/checkdata function to convert the input data to raw
% data (i.e. as if it were coming straight from preprocessing).
% Furthermore, the output data is NOT converted back any more to the
% input data, i.e. the output data is the same as what it would be
% on raw data as input, regardless of the actual input.
%
% Revision 1.38  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.37  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.36  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.35  2007/02/13 16:02:28  roboos
% removed default for cfg.inwardshift, fxed bug in underlying headsurface function, added dipole layer as black points to the 3D figure when cfg.feedback=yes
%
% Revision 1.34  2007/02/13 15:12:51  roboos
% removed cfg.plot3d option
%
% Revision 1.33  2006/10/12 07:23:52  roboos
% changed cfg.headshape='headmodel' into cfg.headshape=[], default behaviour for empty is ok
%
% Revision 1.32  2006/10/09 15:22:28  roboos
% added support for bti148
%
% Revision 1.31  2006/10/03 13:00:50  roboos
% construct dipole grid from inward shifted brain (based on volume model), option 'cortex' is not yet supported by headsurface
%
% Revision 1.30  2006/09/19 16:10:43  roboos
% updated documentation
%
% Revision 1.29  2006/07/24 07:59:16  roboos
% updated documentation
%
% Revision 1.28  2006/07/24 07:51:30  roboos
% improved ducumentation, added two hidden options
% improved memory efficiency in case verify=yes (which is the default)
%
% Revision 1.27  2006/06/14 11:51:06  roboos
% changed a fprintf to be more human-readible
%
% Revision 1.26  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.25  2006/04/10 16:33:46  ingnie
% updated documentation
%
% Revision 1.24  2006/04/06 16:17:10  ingnie
% updated documentation
%
% Revision 1.23  2006/02/24 16:37:56  roboos
% only change in whitespace
%
% Revision 1.22  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.21  2006/01/30 14:34:59  roboos
% added two comments as a reminder that prepare_dipole_grid should be used in the future
%
% Revision 1.20  2006/01/30 14:09:04  roboos
% replaced ctf specific code by call to read_fcdc_header
%
% Revision 1.19  2005/12/14 14:16:38  roboos
% switched from prepare_brain_surface to headsurface subfunction
%
% Revision 1.18  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.17  2005/06/08 13:40:33  roboos
% replaced the specific call to either meg_leadfield or eeg_leadfield to the generic compute_leadfield
%
% Revision 1.16  2005/06/07 13:03:55  roboos
% fixed bug due to change in read_ctf_res4 and ctf2grad
%
% Revision 1.15  2005/06/02 12:18:41  roboos
% changed handling of input data: All input data that contains averages is converted to raw trials (like the output from preprocessing) prior to further processing. The output data is converted back into a format similar to the original input data using RAW2MEG.
%
% Revision 1.14  2005/05/23 09:31:49  roboos
% now gives error if no template is specified
%
% Revision 1.13  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.12  2005/01/17 14:40:38  roboos
% small change in documentation
%
% Revision 1.11  2004/08/20 14:30:27  roboos
% added literatire reference to knosche2002
%
% Revision 1.10  2004/05/05 13:55:20  roberto
% fixed bug that occured when average data was made with keeptrials=yes
%
% Revision 1.9  2004/04/27 13:49:31  roberto
% fixed bug that occurred when template set contained mixture of 150 and 151 channel helmets
%
% Revision 1.8  2004/04/26 09:36:26  roberto
% interp.time was not returned correctly
%
% Revision 1.7  2004/04/22 12:51:19  roberto
% fixed bug that occurred if no rest-channels are present (variable remained undefined)
%
% Revision 1.6  2004/04/13 16:31:09  roberto
% fixed bug in dbstack selection of function filename for Matlab 6.1
%
% Revision 1.5  2004/04/13 14:25:24  roberto
% wrapped code-snippet around mfilename to make it compatible with Matlab 6.1
%
% Revision 1.4  2004/04/08 16:05:19  roberto
% fixed bug in plot3d section where 151 channels were hardcoded
%
% Revision 1.3  2004/04/08 15:45:53  roberto
% added support for datasets with less than 151 MEG channels,
% removed most of the channel selection and now using prepare_vol_sens,
% moved the construction of the brain dipole surface to new function prepare_brain_surface
%

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'ismeg', 'yes');

% set the default configuration 
if ~isfield(cfg, 'headshape'),     cfg.headshape = [];            end
if ~isfield(cfg, 'pruneratio'),    cfg.pruneratio = 1e-3;         end
if ~isfield(cfg, 'spheremesh'),    cfg.spheremesh = 642;          end
if ~isfield(cfg, 'verify'),        cfg.verify = 'yes';            end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'yes';          end
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';            end
if ~isfield(cfg, 'channel'),       cfg.channel = 'MEG';           end
if ~isfield(cfg, 'topoparam'),     cfg.topoparam = 'rms';         end

if ~isfield(cfg, 'inwardshift'),
  % it depends on the volume model and/or headshape that is used for constructing the dipole sheet
  error('you should specify cfg.inwardshift');
end

if isfield(cfg, 'plot3d')
  cfg.feedback = cfg.plot3d;
  cfg = rmfield(cfg, 'plot3d');
  warning('cfg.plot3d is deprecated, see http://www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:development:deprecated for more details.');
end

% for backward compatibility
if ischar(cfg.headshape) && strcmp(cfg.headshape, 'headmodel')
  cfg.headshape = [];
end

if ~isfield(cfg, 'template'), 
  error('you must specify one or more template CTF datasets in cfg.template');
end

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
case 'bti148'
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
orig.dirCF = data.grad.pnt(indxF,:) - data.grad.pnt(indxC,:);
orig.dirRL = data.grad.pnt(indxL,:) - data.grad.pnt(indxR,:);
% construct three orthonormal direction vectors
orig.dirX = normalize(orig.dirCF);
orig.dirY = normalize(orig.dirRL - dot(orig.dirRL, orig.dirX) * orig.dirX);
orig.dirZ = cross(orig.dirX, orig.dirY);
orig.tra = fixedbody(data.grad.pnt(indxC,:), orig.dirX, orig.dirY, orig.dirZ);

% compute the homogenous transformation matrix and transform the positions
tra = inv(templ.tra) * orig.tra;
pnt = data.grad.pnt;
pnt(:,4) = 1;
pnt = (tra * pnt')';
pnt = pnt(:,1:3);

% remove the translation from the transformation matrix and rotate the orientations
tra(:,4) = [0 0 0 1]';
ori = data.grad.ori;
ori(:,4) = 1;
ori = (tra * ori')';
ori = ori(:,1:3);

% construct the final template gradiometer definition
template = [];
template.grad.pnt = pnt;
template.grad.ori = ori;
% keep the same labels and the linear weights of the coils
template.grad.label = data.grad.label;
template.grad.tra   = data.grad.tra;
template.grad.unit  = data.grad.unit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE_VOL_SENS will match the data labels, the gradiometer labels and the 
% volume model labels (in case of a multisphere model) and result in a gradiometer 
% definition that only contains the gradiometers that are present in the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpdat = [];
tmpdat.label = data.label;
tmpcfg = [];
if isfield(cfg, 'hdmfile')
  tmpcfg.hdmfile = cfg.hdmfile;
elseif isfield(cfg, 'vol') 
  tmpcfg.vol = cfg.vol;
end
% note that in the following it is neccessary to keep the two volume
% conduction models seperate, since the single-shell Nolte model contains
% gradiometer specific precomputed parameters
tmpdat.grad  = data.grad;
[volold, data.grad] = prepare_headmodel(tmpcfg, tmpdat);
tmpdat.grad = template.grad;
[volnew, template.grad] = prepare_headmodel(tmpcfg, tmpdat);

% Continue to work with the volume conduction model that matches the
% original data and gradiometer definition. For the template gradiometer
% definition it would not be possible to construct a multi-sphere model.
% For a single sphere model it does not make a difference.

fprintf('mean distance towards template gradiometers is %.2f %s\n', mean(sum((data.grad.pnt-template.grad.pnt).^2, 2).^0.5), template.grad.unit);

% create the dipole grid on which the data will be projected
grid = prepare_dipole_grid(cfg, volold, data.grad);
pos = grid.pos;
clear grid

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
  % also compute the residual variance when interpolating
  rvrealign = rv(data.trial{i}, data.realign{i});
  fprintf('original -> template             RV %.2f %%\n', 100 * mean(rvrealign));
  if strcmp(cfg.verify, 'yes')
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
  tmpcfg.grid.pos = pos;
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
  fprintf('adding %d non-MEG channels back to the data\n', length(rest.label));
  for trial=1:length(rest.trial)
    interp.trial{trial} = [interp.trial{trial}; rest.trial{trial}];
  end
  interp.label = [interp.label; rest.label];
end

% store the configuration of this function call, including that of the previous function call
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: megrealign.m,v 1.49 2008/05/06 15:43:46 sashae Exp $';
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

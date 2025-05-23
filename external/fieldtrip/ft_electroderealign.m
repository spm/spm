function [elec_realigned] = ft_electroderealign(cfg, elec_input)

% FT_ELECTRODEREALIGN rotates, translates, scales, warps and/or projects EEG
% electrode positions. The different methods are described in detail below.
%
% INTERACTIVE - This displays the skin surface together with the electrode or
% gradiometer positions, and allows you to manually adjust (using the graphical user
% interface) the rotation, translation and scaling parameters, so that the electrodes
% correspond with the skin.
%
% FIDUCIAL - This applies a rigid body realignment based on three fiducial or
% anatomical locations. After realigning, the fiducials that are part of the input
% electrode set (typically nose, left and right ear) are along the same axes as
% the fiducials in the target electrode set.
%
% TEMPLATE - This applies a linear or non-linear spatial transformation/deformation
% that automatically minimizes the distance between the input electrodes and a target
% electrode set. The warping methods use a non-linear search to minimize the distance
% between the input electrode positions and the corresponding target electrode.
%
% HEADSHAPE - This applies a spatial transformation/deformation that automatically
% minimizes the distance between the electrodes and the head surface. The warping
% methods use a non-linear search to minimize the distance between the input electrode
% positions and the projection of the electrodes on the head surface.
%
% PROJECT - This projects each of the electrodes to the nearest point on the
% head surface mesh.
%
% MOVEINWARD - This moves all electrodes inward according to their normals.
%
% MNI - This transforms the electrodes nonlinearly using the same transformation of
% the individual anatomical MRI to the MNI template.
%
% Use as
%   [elec_realigned] = ft_electroderealign(cfg)
% with the electrode or gradiometer details in the configuration, or as
%   [elec_realigned] = ft_electroderealign(cfg, elec_input)
% with the electrode or gradiometer definition as 2nd input argument.
%
% The configuration can contain the following options
%   cfg.method         = string representing the method for aligning or placing the electrodes
%                        'interactive'     realign manually using a graphical user interface
%                        'fiducial'        realign using three fiducials (e.g. NAS, LPA and RPA)
%                        'template'        realign the electrodes to match a target electrode set
%                        'headshape'       realign the electrodes to fit the head surface
%                        'project'         projects electrodes onto the head surface
%                        'moveinward'      moves electrodes inward along their normals
%                        'mni'             transforms electrodes from individual subject to MNI space
%   cfg.warp          = string describing the spatial transformation for the template and headshape methods
%                        'rigidbody'       apply a rigid-body warp (default)
%                        'globalrescale'   apply a rigid-body warp with global rescaling
%                        'traditional'     apply a rigid-body warp with individual axes rescaling
%                        'nonlin1'         apply a 1st order non-linear warp
%                        'nonlin2'         apply a 2nd order non-linear warp
%                        'nonlin3'         apply a 3rd order non-linear warp
%                        'nonlin4'         apply a 4th order non-linear warp
%                        'nonlin5'         apply a 5th order non-linear warp
%                        'dykstra2012'     back-project ECoG onto the cortex using energy minimzation
%                        'hermes2010'      back-project ECoG onto the cortex along the local norm vector
%                        'fsaverage'       surface-based realignment with FreeSurfer fsaverage brain (left->left or right->right)
%                        'fsaverage_sym'   surface-based realignment with FreeSurfer fsaverage_sym left hemisphere (left->left or right->left)
%                        'fsinflated'      surface-based realignment with FreeSurfer individual subject inflated brain (left->left or right->right)
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.keepchannel    = string, 'yes' or 'no' (default = 'no')
%   cfg.fiducial       = cell-array with the name of three fiducials used for aligning to the target (default = {'nasion', 'lpa', 'rpa'})
%   cfg.casesensitive  = 'yes' or 'no', determines whether string comparisons between electrode labels are case sensitive (default = 'yes')
%   cfg.feedback       = 'yes' or 'no' (default = 'no')
%
% The electrode positions can be present in the 2nd input argument or can be specified as
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%
% If your input EEG electrodes include the positions of anatomical landmarks or
% fiducials, you can specify the target location of those, for example here in
% millimeter according to the CTF coordinate system with the nose along X and the
% ears along +Y and -Y
%   cfg.target.pos(1,:) = [110 0 0]     % target location of the nose
%   cfg.target.pos(2,:) = [0  90 0]     % target location of the left ear
%   cfg.target.pos(3,:) = [0 -90 0]     % target location of the right ear
%   cfg.target.label    = {'NAS', 'LPA', 'RPA'}
%
% If you want to align the input EEG electrodes to a target electrode sets or to
% multiple target electrode sets (which will be averaged), you should specify the
% target electrode sets either as electrode structures (i.e. when they are already
% read in memory) or as their file names using
%   cfg.target          = single electrode set to serve as the target
% or
%   cfg.target{1..N}    = list of electrode sets to serve as the target, these will be averaged
%
% If you want to align EEG electrodes to the head surface, you should specify the head surface as
%   cfg.headshape      = a filename containing headshape, a structure containing a
%                        single triangulated boundary, or a Nx3 matrix with surface
%                        points
%
% If you want to align ECoG electrodes to the pial surface, you first need to compute
% the cortex hull with FT_PREPARE_MESH. Then use either the algorithm described in
% Dykstra et al. (2012, Neuroimage) or in Hermes et al. (2010, J Neurosci methods) to
% snap the electrodes back to the cortical hull, e.g.
%   cfg.method         = 'headshape'
%   cfg.warp           = 'dykstra2012', or 'hermes2010'
%   cfg.headshape      = a filename containing headshape, a structure containing a
%                        single triangulated boundary, or a Nx3 matrix with surface
%                        points
%   cfg.feedback       = 'yes' or 'no' (default), feedback of the iteration procedure
%
% Additional configuration options for cfg.warp='dykstra2012'
%   cfg.maxiter        = number (default: 50), maximum number of optimization iterations
%   cfg.pairmethod     = 'pos' (default) or 'label', the method for electrode
%                        pairing on which the deformation energy is based
%   cfg.isodistance    = 'yes', 'no' (default) or number, to enforce isotropic
%                        inter-electrode distances (pairmethod 'label' only)
%   cfg.deformweight   = number (default: 1), weight of deformation relative
%                        to shift energy cost (lower increases grid flexibility)
%
% If you want to move the electrodes inward, you should specify
%   cfg.moveinward     = number, the distance that the electrode should be moved
%                        inward (negative numbers result in an outward move)
%
% If you want to align ECoG electrodes to the freesurfer average brain, you should
% specify the path to your headshape (e.g., lh.pial), and ensure you have the
% corresponding registration file (e.g., lh.sphere.reg) in the same directory.
% Moreover, the path to the local freesurfer home is required. Note that, because the
% electrodes are being aligned to the fsaverage brain, the corresponding brain should
% be also used when plotting the data, i.e. use freesurfer/subjects/fsaverage/surf/lh.pial
% rather than surface_pial_left.mat
%   cfg.method         = 'headshape'
%   cfg.warp           = 'fsaverage'
%   cfg.headshape      = string, filename containing subject headshape (e.g. <path to freesurfer/surf/lh.pial>)
%   cfg.fshome         = string, path to freesurfer
%
% If you want to transform electrodes from individual subject coordinates to MNI
% space, you should specify the following
%   cfg.mri            = structure with the individual anatomical MRI relative to which electrodes are specified, or the filename of the MRI, see FT_READ_MRI
%   cfg.templatemri    = string, filename of the MNI template (default = 'T1.mnc' for SPM2 or 'T1.nii' for SPM8 and SPM12)
%   cfg.spmversion     = string, 'spm2', 'spm8', 'spm12' (default = 'spm12')
%   cfg.spmmethod      = string, 'old', 'new' or 'mars', see FT_VOLUMENORMALISE
%   cfg.nonlinear      = string, 'yes' or 'no', see FT_VOLUMENORMALISE
%
% See also FT_READ_SENS, FT_VOLUMEREALIGN, FT_INTERACTIVEREALIGN,
% FT_DETERMINE_COORDSYS, FT_PREPARE_MESH

% Copyright (C) 2005-2025, Robert Oostenveld, Arjen Stolk
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% the interactive method uses a global variable to get the data from the figure when it is closed
global realigned

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    elec_input
ft_preamble provenance elec_input

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'forbidden',  {'outline'});
cfg = ft_checkconfig(cfg, 'renamed',    {'template', 'target'});
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'realignfiducials', 'fiducial'});
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'realignfiducial',  'fiducial'});
cfg = ft_checkconfig(cfg, 'renamedval', {'warp', 'homogenous', 'rigidbody'});
cfg = ft_checkconfig(cfg, 'renamedval', {'warp', 'homogeneous', 'rigidbody'});

% set these defaults already here, the rest will follow later
cfg.target        = ft_getopt(cfg, 'target',  []);       % for electrodes or fiducials, always with labels
cfg.coordsys      = ft_getopt(cfg, 'coordsys');          % this allows for automatic template fiducial placement
if ~isempty(cfg.coordsys) && isempty(cfg.target)
  % set the template fiducial locations according to the coordinate system
  switch lower(cfg.coordsys)
    case 'ctf'
      cfg.target = [];
      cfg.target.coordsys = 'ctf';
      cfg.target.pos(1,:) = [100  0 0];
      cfg.target.pos(2,:) = [0   80 0];
      cfg.target.pos(3,:) = [0  -80 0];
      cfg.target.label{1} = 'NAS';
      cfg.target.label{2} = 'LPA';
      cfg.target.label{3} = 'RPA';
    otherwise
      ft_error('the %s coordinate system is not automatically supported, please specify fiducial details in cfg.target')
  end
elseif isempty(cfg.target)
  % remove the field, otherwise ft_checkconfig will complain depending on
  % the settings for checkconfig
  cfg = rmfield(cfg, 'target');
end

% ensure that the right cfg options have been set corresponding to the method
switch cfg.method
  case 'template'      % realign the sensors to match a template set
    cfg = ft_checkconfig(cfg, 'required', 'target', 'forbidden', 'headshape');
  case 'headshape'     % realign the sensors to fit the head surface
    cfg = ft_checkconfig(cfg, 'required', 'headshape', 'forbidden', 'target');
  case 'fiducial'      % realign using the NAS, LPA and RPA fiducials
    cfg = ft_checkconfig(cfg, 'required', 'target', 'forbidden', 'headshape');
  case 'moveinward'    % moves eletrodes inward
    cfg = ft_checkconfig(cfg, 'required', 'moveinward');
end % switch cfg.method

% set the rest of the defaults, this is to avoid confusion in ft_checkconfig w.r.t. to the method specific checks
cfg.warp          = ft_getopt(cfg, 'warp', 'rigidbody');
cfg.channel       = ft_getopt(cfg, 'channel',  'all');
cfg.keepchannel   = ft_getopt(cfg, 'keepchannel', 'no');
cfg.feedback      = ft_getopt(cfg, 'feedback', 'no');
cfg.casesensitive = ft_getopt(cfg, 'casesensitive', 'no');
cfg.headshape     = ft_getopt(cfg, 'headshape', []);     % for triangulated head surface, without labels

if strcmp(cfg.method, 'fiducial') && isfield(cfg, 'warp') && ~isequal(cfg.warp, 'rigidbody')
  ft_warning('The method ''fiducial'' implies a rigid body tramsformation. See also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1722');
  cfg.warp = 'rigidbody';
end

% the data can be passed as input arguments or can be read from disk
hasdata = exist('elec_input', 'var');

% get the electrode definition that should be warped
if ~hasdata
  elec_input = ft_fetch_sens(cfg);
else
  % the input electrodes were specified as second input argument
  % or read from cfg.inputfile
end

% ensure that the units are specified
elec_input = ft_determine_units(elec_input);

% ensure up-to-date sensor description (Oct 2011)
elec_input = ft_datatype_sens(elec_input);

% ensure that channel and electrode positions are the same
assert(isequaln(elec_input.elecpos, elec_input.chanpos), 'this function requires same electrode and channel positions');

% remember the original electrode locations and labels and do all the work with a
% temporary copy, this involves channel selection and changing to lower case
elec = elec_input;

% instead of working with all sensors, only work with the fiducials
% this is useful for gradiometer structures
if strcmp(cfg.method, 'fiducial') && isfield(elec, 'fid')
  fprintf('using the fiducials instead of the electrode positions\n');
  elec.fid.unit = elec.unit; % copy the units over
  elec          = elec.fid;  
end

usetarget    = isfield(cfg, 'target')    && ~isempty(cfg.target);
useheadshape = isfield(cfg, 'headshape') && ~isempty(cfg.headshape);

if usetarget
  % get the template electrode definitions
  if ~iscell(cfg.target)
    cfg.target = {cfg.target};
  end
  Ntarget = length(cfg.target);
  for i=1:Ntarget
    if isstruct(cfg.target{i})
      target(i) = cfg.target{i};
    else
      target(i) = ft_read_sens(cfg.target{i}, 'senstype', 'eeg');
    end
  end
  clear tmp
  for i=1:Ntarget
    try
      % ensure up-to-date sensor description and that the units are consistent with the electrodes
      tmp(i) = ft_convert_units(ft_datatype_sens(target(i)), elec.unit);
    catch
      ft_warning('cannot check the consistency of the units in the fiducials and the electrodes, you need to ensure this yourself');
      tmp(i) = target(i);
    end
  end
  target = tmp;
end

if useheadshape
  % get the surface describing the head shape
  [headshape.pos, headshape.tri] = headsurface([], [], 'headshape', cfg.headshape);
  if isfield(cfg.headshape, 'unit')
    headshape.unit = cfg.headshape.unit;
  end
  
  % ensure that the units are consistent with the electrodes
  headshape = ft_convert_units(headshape, elec.unit);
end

% convert all labels to lower case for string comparisons
cfg.channel = ft_channelselection(cfg.channel, elec.label);
if strcmp(cfg.casesensitive, 'no')
  elec.label  = lower(elec.label);
  cfg.channel = lower(cfg.channel);
  if usetarget
    for j=1:length(target)
      for i=1:length(target(j).label)
        target(j).label{i} = lower(target(j).label{i});
      end
    end
  end
end
[cfgsel, datsel] = match_str(cfg.channel, elec.label);
% keep the original channel labels
label_original = elec_input.label(datsel);

% start with an empty structure, this will be returned at the end
realigned = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(cfg.method, 'template')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % determine electrode selection and overlapping subset for warping
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  for i=1:Ntarget
    cfg.channel = ft_channelselection(cfg.channel, target(i).label);
  end
  
  % make consistent subselection of electrodes
  [cfgsel, datsel] = match_str(cfg.channel, elec.label);
  elec.label = elec.label(datsel);
  elec.elecpos   = elec.elecpos(datsel,:);
  for i=1:Ntarget
    [cfgsel, datsel] = match_str(cfg.channel, target(i).label);
    target(i).label   = target(i).label(datsel);
    target(i).elecpos = target(i).elecpos(datsel,:);
  end
  
  % compute the average of the target electrode positions
  average = ft_average_sens(target);
  
  fprintf('warping electrodes to average template... '); % the newline comes later
  [realigned.elecpos, realigned.m] = ft_warp_optim(elec.elecpos, average.elecpos, cfg.warp);
  realigned.label = elec.label;
  
  dpre  = mean(sqrt(sum((average.elecpos - elec.elecpos).^2, 2)));
  dpost = mean(sqrt(sum((average.elecpos - realigned.elecpos).^2, 2)));
  fprintf('mean distance prior to warping %f, after warping %f\n', dpre, dpost);
  
  if strcmp(cfg.feedback, 'yes')
    % create an empty figure, continued below...
    figure
    axis equal
    axis vis3d
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    % plot all electrodes before warping
    ft_plot_sens(elec, 'r*');
    
    % plot all electrodes after warping
    ft_plot_sens(realigned, 'm.', 'label', 'label');
    
    % plot the template electrode locations
    ft_plot_sens(average, 'b.');
    
    % plot lines connecting the input electrode locations with the template locations
    % plot lines connecting the realigned electrode locations with the template locations
    my_line3(elec.elecpos, average.elecpos, 'color', 'r');
    my_line3(realigned.elecpos, average.elecpos, 'color', 'm');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'headshape')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % determine electrode selection and overlapping subset for warping
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  [cfgsel, datsel] = match_str(cfg.channel, elec.label);
  elec.label   = elec.label(datsel);
  elec.elecpos = elec.elecpos(datsel,:);
  
  realigned.label = elec.label;
  if strcmp(cfg.warp, 'dykstra2012')
    realigned.elecpos = warp_dykstra2012(cfg, elec, headshape);
  elseif strcmp(cfg.warp, 'hermes2010')
    realigned.elecpos = warp_hermes2010(cfg, elec, headshape);
  elseif strcmp(cfg.warp, 'fsaverage')
    realigned.elecpos = warp_fsaverage(cfg, elec);
  elseif strcmp(cfg.warp, 'fsaverage_sym')
    realigned.elecpos = warp_fsaverage_sym(cfg, elec);
  elseif strcmp(cfg.warp, 'fsinflated')
    realigned.elecpos = warp_fsinflated(cfg, elec);
  else
    fprintf('warping electrodes to skin surface... '); % the newline comes later
    [realigned.elecpos, realigned.m] = ft_warp_optim(elec.elecpos, headshape, cfg.warp);
    
    dpre  = ft_warp_error([],     elec.elecpos, headshape, cfg.warp);
    dpost = ft_warp_error(realigned.m, elec.elecpos, headshape, cfg.warp);
    fprintf('mean distance prior to warping %f, after warping %f\n', dpre, dpost);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'fiducial')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % the fiducials have to be present in the electrodes and in the template set
  label = intersect(lower(elec.label), lower(target.label));
  
  if ~isfield(cfg, 'fiducial') || isempty(cfg.fiducial)
    % try to determine the names of the fiducials automatically
    option1 = {'nasion' 'left' 'right'};
    option2 = {'nasion' 'lpa' 'rpa'};
    option3 = {'nz' 'left' 'right'};
    option4 = {'nz' 'lpa' 'rpa'};
    option5 = {'nas' 'left' 'right'};
    option6 = {'nas' 'lpa' 'rpa'};
    if length(match_str(label, option1))==3
      cfg.fiducial = option1;
    elseif length(match_str(label, option2))==3
      cfg.fiducial = option2;
    elseif length(match_str(label, option3))==3
      cfg.fiducial = option3;
    elseif length(match_str(label, option4))==3
      cfg.fiducial = option4;
    elseif length(match_str(label, option5))==3
      cfg.fiducial = option5;
    elseif length(match_str(label, option6))==3
      cfg.fiducial = option6;
    else
      ft_error('could not determine consistent fiducials in the input and the target, please specify cfg.fiducial or cfg.coordsys')
    end
  end
  fprintf('matching fiducials {''%s'', ''%s'', ''%s''}\n', cfg.fiducial{1}, cfg.fiducial{2}, cfg.fiducial{3});
  
  % determine electrode selection
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  [cfgsel, datsel] = match_str(cfg.channel, elec.label);
  elec.label     = elec.label(datsel);
  elec.elecpos   = elec.elecpos(datsel,:);
  
  if length(cfg.fiducial)~=3
    ft_error('you must specify exactly three fiducials');
  end
  
  % do case-insensitive search for fiducial locations
  nas_indx = match_str(lower(elec.label), lower(cfg.fiducial{1}));
  lpa_indx = match_str(lower(elec.label), lower(cfg.fiducial{2}));
  rpa_indx = match_str(lower(elec.label), lower(cfg.fiducial{3}));
  if length(nas_indx)~=1 || length(lpa_indx)~=1 || length(rpa_indx)~=1
    ft_error('not all fiducials were found in the electrode set');
  end
  elec_nas = elec.elecpos(nas_indx,:);
  elec_lpa = elec.elecpos(lpa_indx,:);
  elec_rpa = elec.elecpos(rpa_indx,:);
  
  % FIXME change the flow in the remainder
  % if one or more template electrode sets are specified, then align to the average of those
  % if no template is specified, then align so that the fiducials are along the axis
  
  % find the matching fiducials in the template and average them
  tmpl_nas = nan(Ntarget,3);
  tmpl_lpa = nan(Ntarget,3);
  tmpl_rpa = nan(Ntarget,3);
  for i=1:Ntarget
    nas_indx = match_str(lower(target(i).label), lower(cfg.fiducial{1}));
    lpa_indx = match_str(lower(target(i).label), lower(cfg.fiducial{2}));
    rpa_indx = match_str(lower(target(i).label), lower(cfg.fiducial{3}));
    if length(nas_indx)~=1 || length(lpa_indx)~=1 || length(rpa_indx)~=1
      ft_error('not all fiducials were found in template %d', i);
    end
    tmpl_nas(i,:) = target(i).elecpos(nas_indx,:);
    tmpl_lpa(i,:) = target(i).elecpos(lpa_indx,:);
    tmpl_rpa(i,:) = target(i).elecpos(rpa_indx,:);
  end
  tmpl_nas = mean(tmpl_nas,1);
  tmpl_lpa = mean(tmpl_lpa,1);
  tmpl_rpa = mean(tmpl_rpa,1);
  
  % realign both to a common coordinate system
  elec2common  = ft_headcoordinates(elec_nas, elec_lpa, elec_rpa);
  templ2common = ft_headcoordinates(tmpl_nas, tmpl_lpa, tmpl_rpa);
  
  % compute the combined transform
  realigned         = [];
  realigned.m       = templ2common \ elec2common;
  
  % apply the transformation to the fiducials as sanity check
  realigned.elecpos(1,:) = ft_warp_apply(realigned.m, elec_nas, 'homogeneous');
  realigned.elecpos(2,:) = ft_warp_apply(realigned.m, elec_lpa, 'homogeneous');
  realigned.elecpos(3,:) = ft_warp_apply(realigned.m, elec_rpa, 'homogeneous');
  realigned.label        = cfg.fiducial;
  
  nas_indx = match_str(lower(elec.label), lower(cfg.fiducial{1}));
  lpa_indx = match_str(lower(elec.label), lower(cfg.fiducial{2}));
  rpa_indx = match_str(lower(elec.label), lower(cfg.fiducial{3}));
  dpre  = mean(sqrt(sum((elec.elecpos([nas_indx lpa_indx rpa_indx],:) - [tmpl_nas; tmpl_lpa; tmpl_rpa]).^2, 2)));
  nas_indx = match_str(lower(realigned.label), lower(cfg.fiducial{1}));
  lpa_indx = match_str(lower(realigned.label), lower(cfg.fiducial{2}));
  rpa_indx = match_str(lower(realigned.label), lower(cfg.fiducial{3}));
  dpost = mean(sqrt(sum((realigned.elecpos([nas_indx lpa_indx rpa_indx],:) - [tmpl_nas; tmpl_lpa; tmpl_rpa]).^2, 2)));
  fprintf('mean distance between fiducials prior to realignment %f, after realignment %f\n', dpre, dpost);
  
  if strcmp(cfg.feedback, 'yes')
    % create an empty figure, continued below...
    figure
    axis equal
    axis vis3d
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    % plot the first three electrodes before transformation
    my_plot3(elec.elecpos(1,:), 'r*');
    my_plot3(elec.elecpos(2,:), 'r*');
    my_plot3(elec.elecpos(3,:), 'r*');
    my_text3(elec.elecpos(1,:), elec.label{1}, 'color', 'r');
    my_text3(elec.elecpos(2,:), elec.label{2}, 'color', 'r');
    my_text3(elec.elecpos(3,:), elec.label{3}, 'color', 'r');
    
    % plot the template fiducials
    my_plot3(tmpl_nas, 'b*');
    my_plot3(tmpl_lpa, 'b*');
    my_plot3(tmpl_rpa, 'b*');
    my_text3(tmpl_nas, ' nas', 'color', 'b');
    my_text3(tmpl_lpa, ' lpa', 'color', 'b');
    my_text3(tmpl_rpa, ' rpa', 'color', 'b');
    
    % plot all electrodes after transformation
    my_plot3(realigned.elecpos, 'm.');
    my_plot3(realigned.elecpos(1,:), 'm*');
    my_plot3(realigned.elecpos(2,:), 'm*');
    my_plot3(realigned.elecpos(3,:), 'm*');
    my_text3(realigned.elecpos(1,:), realigned.label{1}, 'color', 'm');
    my_text3(realigned.elecpos(2,:), realigned.label{2}, 'color', 'm');
    my_text3(realigned.elecpos(3,:), realigned.label{3}, 'color', 'm');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'interactive')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  tmpcfg = [];
  tmpcfg.individual.elec = elec;
  if useheadshape
    tmpcfg.template.headshape = headshape;
  end
  if usetarget
    if iscell(target)
      if numel(target)>1
        ft_notice('computing the average electrode positions');
        tmpcfg.template.elec = ft_average_sens(target);
      else
        tmpcfg.template.elec = target{1};
      end
    elseif isstruct(cfg.target)
      tmpcfg.template.elec = target;
    end
    tmpcfg.template.elecstyle = {'facecolor', 'blue'};
    ft_info('plotting the target electrodes in blue');
  end
  
  % use the more generic ft_interactiverealign for the actual work
  tmpcfg = ft_interactiverealign(tmpcfg);
  % only keep the transformation, it will be applied to the electrodes further down
  realigned.m = tmpcfg.m;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'project')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % determine electrode selection
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  [cfgsel, datsel] = match_str(cfg.channel, elec.label);
  elec.label     = elec.label(datsel);
  elec.elecpos   = elec.elecpos(datsel,:);
  
  realigned.label = elec.label;
  [dum, realigned.elecpos] = project_elec(elec.elecpos, headshape.pos, headshape.tri);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'moveinward')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % determine electrode selection
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  [cfgsel, datsel] = match_str(cfg.channel, elec.label);
  elec.label     = elec.label(datsel);
  elec.elecpos   = elec.elecpos(datsel,:);
  
  realigned.label = elec.label;
  realigned.elecpos = surface_shift(elec.elecpos, [], -cfg.moveinward); % move inward with the specified amount, hence negative
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'mni')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % get the anatomical mri
  if ischar(cfg.mri)
    mri = ft_read_mri(cfg.mri);
  else
    mri = cfg.mri;
  end
  
  % ensure consistent units of the electrodes and individual anatomical MRI with the template MRI (assumed to be in mm)
  mri  = ft_convert_units(mri,  'mm');
  elec = ft_convert_units(elec, 'mm');
  
  % spatial normalisation of the MRI to the template
  tmpcfg           = keepfields(cfg, {'spmversion', 'spmmethod', 'nonlinear'});
  if isfield(cfg, 'templatemri')
    % this option is called differently for the two functions
    tmpcfg.template = cfg.templatemri;
  end
  normalise = ft_volumenormalise(tmpcfg, mri);
  
  % the normalisation from original subject head coordinates to MNI consists of an initial rigid body transformation, followed by a more precise (non)linear transformation
  realigned.label = elec.label;
  realigned.elecpos = ft_warp_apply(normalise.params, ft_warp_apply(normalise.initial, elec.elecpos, 'homogeneous'), 'individual2sn');
  
else
  ft_error('unknown method');
end % if method


% apply the spatial transformation to the original input electrodes, and replace the
% electrode labels by their case-sensitive original labels
switch cfg.method
  case {'template', 'headshape'}
    if strcmpi(cfg.warp, 'dykstra2012') || strcmpi(cfg.warp, 'hermes2010') || ...
        strcmpi(cfg.warp, 'fsaverage') || strcmpi(cfg.warp, 'fsaverage_sym') || strcmpi(cfg.warp, 'fsinflated')
      elec_realigned = realigned;
      elec_realigned.label = label_original;
    else
      % the transformation is a linear or non-linear warp, i.e. a vector
      try
        % convert the vector with fitted parameters into a 4x4 homogenous transformation
        % apply the transformation to the original complete set of sensors
        elec_realigned = ft_transform_geometry(feval(cfg.warp, realigned.m), elec_input);
      catch
        % the previous section will fail for nonlinear transformations
        elec_realigned.label = label_original;
        try
          elec_realigned.elecpos = ft_warp_apply(realigned.m, elec_input.elecpos, cfg.warp);
        end % FIXME why is an error here not dealt with?
      end
      % remember the transformation
      elec_realigned.(cfg.warp) = realigned.m;
    end
    
  case  {'fiducial', 'interactive'}
    % the transformation is a 4x4 homogenous matrix
    % apply the transformation to the original complete set of sensors
    elec_realigned = ft_transform_geometry(realigned.m, elec_input);
    % remember the transformation
    elec_realigned.homogeneous = realigned.m;
    
  case {'project', 'moveinward', 'mni'}
    % nothing to be done
    elec_realigned = realigned;
    elec_realigned.label = label_original;
    
  otherwise
    ft_error('unknown method');
end % switch method

% the coordinate system is in general not defined after transformation
if isfield(elec_realigned, 'coordsys')
  elec_realigned = rmfield(elec_realigned, 'coordsys');
end

% in some cases the coordinate system matches that of the input target or headshape
switch cfg.method
  case 'template'
    if isfield(target, 'coordsys')
      elec_realigned.coordsys = target.coordsys;
    end
  case 'headshape'
    if isfield(headshape, 'coordsys')
      elec_realigned.coordsys = headshape.coordsys;
    end
    if isfield(elec_input, 'coordsys')
      if strcmp(cfg.warp, 'dykstra2012') || strcmp(cfg.warp, 'hermes2010')  % this warp simply moves the electrodes in the same coordinate space
        elec_realigned.coordsys = elec_input.coordsys;
      elseif strcmp(cfg.warp, 'fsaverage')
        elec_realigned.coordsys = 'fsaverage';
      elseif strcmp(cfg.warp, 'fsaverage_sym')
        elec_realigned.coordsys = 'fsaverage_sym';
      end
    end
  case 'fiducial'
    if isfield(target, 'coordsys')
      elec_realigned.coordsys = target.coordsys;
    end
  case 'interactive'
    % the coordinate system is not known
  case {'project', 'moveinward'}
    % the coordinate system remains the same
    if isfield(elec_input, 'coordsys')
      elec_realigned.coordsys = elec_input.coordsys;
    end
  case {'mni'}
    % the coordinate system remains the same
    if isfield(normalise, 'coordsys')
      elec_realigned.coordsys = normalise.coordsys;
    else
      elec_realigned.coordsys = 'mni';
    end
  otherwise
    ft_error('unknown method');
end

if istrue(cfg.keepchannel)
  % append the channels that are not realigned
  [dum, idx] = setdiff(elec_input.label, elec_realigned.label);
  idx = sort(idx);
  elec_realigned.label = [elec_realigned.label; elec_input.label(idx)];
  elec_realigned.elecpos = [elec_realigned.elecpos; elec_input.elecpos(idx,:)];
end

% copy over unit, chantype, chanunit, and tra information in case this was not already done
if ~isfield(elec_realigned, 'unit') && isfield(elec_input, 'unit')
  elec_realigned.unit = elec_input.unit;
end
if ~isfield(elec_realigned, 'chantype') && isfield(elec_input, 'chantype')
  idx = match_str(elec_input.label, elec_realigned.label);
  elec_realigned.chantype = elec_input.chantype(idx);
end
if ~isfield(elec_realigned, 'chanunit') && isfield(elec_input, 'chanunit')
  elec_realigned.chanunit = elec_input.chanunit;
  idx = match_str(elec_input.label, elec_realigned.label);
  elec_realigned.chanunit = elec_input.chanunit(idx);
end

% update it to the latest version
elec_realigned = ft_datatype_sens(elec_realigned);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   elec_input
ft_postamble provenance elec_realigned
ft_postamble history    elec_realigned

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some simple SUBFUNCTIONs that facilitate 3D plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = my_plot3(xyz, varargin)
h = plot3(xyz(:,1), xyz(:,2), xyz(:,3), varargin{:});
function h = my_text3(xyz, varargin)
h = text(xyz(:,1), xyz(:,2), xyz(:,3), varargin{:});
function my_line3(xyzB, xyzE, varargin)
for i=1:size(xyzB,1)
  line([xyzB(i,1) xyzE(i,1)], [xyzB(i,2) xyzE(i,2)], [xyzB(i,3) xyzE(i,3)], varargin{:})
end

function [vol, sens, cfg] = prepare_headmodel(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that helps to prepare the electrodes/gradiometers and the volume
% this is used in sourceanalysis and dipolefitting
%
% This function will get the gradiometer/electrode definition from 
%   cfg.channel
%   cfg.gradfile
%   cfg.elecfile
%   cfg.elec
%   cfg.grad
%   data.grad
%   data.elec
% and the volume conductor definition from 
%   cfg.hdmfile
%   cfg.vol
%   data.vol
% Subsequently it will remove the gradiometers/electrodes that are not
% present in the data. Finally it with attach the gradiometers to a
% multi-sphere head model (if supplied) or attach the electrodes to
% a BEM head model.
% 
% This function will return the electrodes/gradiometers in an order that is
% consistent with the original order in the input electrode/gradiometer file
% or definition. This is also reflected in sens.label.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2004-2008, Robert Oostenveld
%
% $Log: prepare_headmodel.m,v $
% Revision 1.2  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.1  2008/04/10 07:57:37  roboos
% renamed from prepare_vol_sens into prepare_headmodel, based on rev
% 1.28. This is to avild a name clash with the lower level function
% that is part of the forwinv module (i.e. this is part of the
% rearrangement of high-level and low-level functionality between
% fieldtrip and teh seperate modules)
%
% Revision 1.28  2008/03/18 12:26:47  roboos
% use senstype function to add descriptive string to the output (sens.type=...)
%
% Revision 1.27  2008/03/05 10:50:06  roboos
% switched to read_sens for reading elec or grad structures from file
% moved convert_ama2vol to standalone function
%
% Revision 1.26  2007/04/19 17:15:56  roboos
% only initialize the nolte method when the gradiometer array is not empty (usefull for plotting)
%
% Revision 1.25  2006/07/20 15:04:23  roboos
% added instructive fprintf message
%
% Revision 1.24  2006/05/23 10:17:32  roboos
% changed some comments, no code changes
%
% Revision 1.23  2006/04/12 08:38:15  ingnie
% updated documentation
%
% Revision 1.22  2006/04/10 16:35:20  ingnie
% updated documentation
%
% Revision 1.21  2006/03/21 09:44:03  roboos
% implemented support for Guido Nolte's method using vol.type='nolte'
% for neuromag: store the surface normals in the vol.bnd.nrm (consistent with nolte)
%
% Revision 1.20  2005/12/14 10:42:26  roboos
% removed warning for synthetic gradiometers
% if topolabel is present in data, look at that instead of label (for ICA)
% always try to add vol.skin and vol.brain
%
% Revision 1.19  2005/11/16 09:14:06  roboos
% moved the conversion from dipoli/ama format to fieldtrip/vol format to subfunction
% changed the post-processing of the EEG BEM model: incorporate the tra and mat matrices into one matrix (for computational efficiency)
%
% Revision 1.18  2005/09/29 01:15:32  roboos
% changed construction of the localsphere per coil for CTF hdm file
%
% Revision 1.17  2005/09/29 00:47:05  roboos
% added support for mbfys_ama BEM volume conductor file
%
% Revision 1.16  2005/08/05 07:26:22  roboos
% switched to teh use of read_fcdc_elec for cfg.gradfile
%
% Revision 1.15  2005/07/18 10:13:23  roboos
% fixed reading gradiometer from cfg.gradfile for CTF using ctf2grad
% added reading gradiometer from fif file
%
% Revision 1.14  2005/06/28 19:50:42  roboos
% minor change to solve problem when compiling to c-code, functionality is the same
% load(cfg.filename) gives compilation error, so first copy the string with the filename from the structure in a plain variable and then load(...)
%
% Revision 1.13  2005/06/08 16:33:54  roboos
% changed reading of gradiometers from matlab file
% prune the coils for a higher-order gradiometer system (i.e. remove non-contributing coils) -> this temporary solves a problem with multisphere headmodels
%
% Revision 1.12  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.11  2005/03/03 10:52:37  roboos
% changed the handling of the channel selection. The input cfg now specifies
% the channels that should be kept in the sensor array. If no selection is
% specified in the input cfg (e.g. for dipolesimulation), the default is 'all'.
%
% Revision 1.10  2005/01/26 08:05:23  roboos
% added empty data if second input argument not given
%
% Revision 1.9  2005/01/17 14:52:50  roboos
% changed to use read_fcdc_elec
%
% Revision 1.8  2004/12/08 18:00:13  roboos
% implemented consistent method of selecting a subset of channels for
% forward and inverse computations using cfg.channel and updated the
% ducumentation
%
% Revision 1.7  2004/09/06 08:45:38  roboos
% moved reading of neuromag BEM bondary from find_inside_vol into prepare_vol_sens
%
% Revision 1.6  2004/09/03 09:15:28  roboos
% added channelselection to the volume for neuromag
%
% Revision 1.5  2004/09/03 06:39:09  roboos
% fixed bug in non-functional neuromag section, cfg->vol
%
% Revision 1.4  2004/09/01 17:30:09  roboos
% added some explanation to the chansel for neuromag
%
% Revision 1.3  2004/08/31 13:55:22  roboos
% added initialization of neuromag megmodel
%
% Revision 1.2  2004/08/06 08:54:53  roboos
% fixed bug for gradiometer info coming from matlab file
%
% Revision 1.1  2004/06/28 08:59:38  roboos
% moved files from fieldtrip to fieldtrip/private
%
% Revision 1.2  2004/05/08 21:06:25  roberto
% added support for reading gradiometer from matlab file
%
% Revision 1.1  2004/04/08 15:51:20  roberto
% initial submissiion into cvs
%

% set the defaults
if ~isfield(cfg, 'channel'), cfg.channel = 'all';   end
if ~isfield(cfg, 'order'),   cfg.order = 10;        end % order of expansion for Nolte method; 10 should be enough for real applications; in simulations it makes sense to go higher

if nargin<2
  data = [];
elseif isfield(data, 'topolabel')
  % the data reflects a componentanalysis, where the topographic and the
  % timecourse labels are different
  cfg.channel = channelselection(cfg.channel, data.topolabel);
elseif isfield(data, 'label')
  % In the subsequent code, the matching channels in the sensor array and
  % in the configuration will be selected. To ensure that these channels
  % are also present in the data, update the configuration to match the data.
  cfg.channel = channelselection(cfg.channel, data.label);
end

% get the volume conduction model
if isfield(cfg, 'hdmfile')
  fprintf('reading headmodel from file %s\n', cfg.hdmfile);
  switch(filetype(cfg.hdmfile))
    case 'matlab'
      matfile = cfg.hdmfile;   % this solves a problem with the matlab compiler v3
      vol = getfield(load(matfile), 'vol');
    
    case 'ctf_hdm'
      vol = read_ctf_hdm(cfg.hdmfile);
    
    case 'asa_vol'
      vol = read_asa_vol(cfg.hdmfile);
    
    case 'mbfys_ama'
      ama = loadama(cfg.hdmfile);
      vol = ama2vol(ama);
  
    case 'neuromag_fif'
      % do not read the volume into Matlab, but use external Neuromag toolbox
      vol.type     = 'neuromag';
      vol.filename = cfg.hdmfile;
      vol.chansel  = [];  % this is defined later based on the channels present in the data
      % initialize the Neuromag toolbox, this requires a gradfile and hdmfile
      fprintf('using Neuromag volume conductor from %s\n', cfg.hdmfile);
      fprintf('using Neuromag gradiometer definition from %s\n', cfg.gradfile);
      megmodel('head', cfg.gradfile, cfg.hdmfile);
      % read the triangulated boundary from the neuromag BEM model
      [vol.bnd.pnt, vol.bnd.tri, vol.bnd.nrm] = loadtri(vol.filename);
      vol.bnd.pnt = vol.bnd.pnt*100;  % convert to cm
    otherwise
      error('unknown fileformat for volume conductor model');
  end
elseif isfield(cfg, 'vol')
  fprintf('using headmodel specified in the configuration\n');
  vol = cfg.vol;
elseif isfield(data, 'vol') 
  fprintf('using headmodel specified in the data\n');
  vol = data.vol;
else
  error('no headmodel specified');
end

% determine the skin compartment
if ~isfield(vol, 'skin')
  if isfield(vol, 'bnd')
    vol.skin   = find_outermost_boundary(vol.bnd);
  elseif isfield(vol, 'r') && length(vol.r<=4)
    [dum, vol.skin] = max(vol.r);
  end
end

% determine the brain compartment
if ~isfield(vol, 'brain')
  if isfield(vol, 'bnd')
    vol.brain  = find_innermost_boundary(vol.bnd);
  elseif isfield(vol, 'r') && length(vol.r<=4)
    [dum, vol.brain] = min(vol.r);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the gradiometer positions and orientations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg, 'gradfile')
  fprintf('reading gradiometers from file %s\n', cfg.gradfile);
  grad = read_sens(cfg.gradfile);
elseif isfield(cfg, 'grad') 
  fprintf('using gradiometers specified in the configuration\n');
  grad = cfg.grad;
elseif isfield(data, 'grad') 
  fprintf('using gradiometers specified in the data\n');
  grad = data.grad;
end

% do postprocessing of the volume and gradiometers in case of MEG
if exist('grad', 'var')
  % update the channelselection with those present in the gradiometer array
  cfg.channel = channelselection(cfg.channel, grad.label);
  % test whether it is a magnetometer or a first-order CTF gradiometer 
  if length(grad.label)==size(grad.pnt,1)
    sel = match_str(grad.label, cfg.channel);
    fprintf('selecting %d magnetometers\n', length(sel));
    grad.label = grad.label(sel);
    grad.pnt   = grad.pnt(sel,:);
    grad.ori   = grad.ori(sel,:);
    if isfield(grad, 'tra')
      grad.tra = sparse(grad.tra(sel, sel));
    end
  elseif length(grad.label)==size(grad.pnt,1)/2
    sel = match_str(grad.label, cfg.channel);
    fprintf('selecting %d first-order gradiometers\n', length(sel));
    % extend the selection to both bottom and top coil
    grad.label = grad.label(sel);
    sel2 = [sel; sel+size(grad.pnt,1)/2];
    grad.pnt   = grad.pnt(sel2,:);
    grad.ori   = grad.ori(sel2,:);
    if isfield(grad, 'tra')
      grad.tra = sparse(grad.tra(sel, sel2));
    end
  else
    sel = match_str(grad.label, cfg.channel);
    grad.label = grad.label(sel);
    % In this case it is not possible to reduce the number of coils, i.e.
    % to make a selection out of pnt and ori. Therefore keep them all and
    % only make a selection out of the linear transformation matrix tra.
    if isfield(grad, 'tra')
      grad.tra = sparse(grad.tra(sel, :));
    else
      error('no transfer matrix present in the gradiometer definition, cannot make channel selection');
    end
    % remove the coils from the grad.pnt and ori field that do not contribute to any channel's output
    selcoil = find(sum(grad.tra,1)~=0);
    grad.pnt = grad.pnt(selcoil,:);
    grad.ori = grad.ori(selcoil,:);
    grad.tra = grad.tra(:,selcoil);
  end
  % If the volume conduction model consists of multiple spheres then we
  % have to match the channels in the gradiometer array and the volume
  % conduction model.
  if isfield(vol, 'label') && length(vol.label)>1
    % get the local spheres for the MEG channels, this will be ordered
    % according to the ordering of the gradiometer channels
    [selgrad, selvol] = match_str(grad.label, vol.label);
    % the CTF way of storing the headmodel is one-sphere-per-channel
    % whereas the FieldTrip way is one-sphere-per-coil  
    Nchans = size(grad.tra,1);
    Ncoils = size(grad.tra,2);
    % for each coil in the MEG helmet, determine the corresponding local sphere
    for i=1:Ncoils
      coilindex = find(grad.tra(:,i)~=0); % to which channel does the coil belong
      if length(coilindex)>1
        % this indicates that there are multiple channels to which this coil contributes
        % which means that grad.tra describes a synthetic higher-order gradient
        error('synthetic gradients not supported during volume conductor setup');
      end
      coillabel = grad.label{coilindex};  % what is the label of the channel
      chanindex = strmatch(coillabel, vol.label, 'exact');
      multisphere.r(i,:) = vol.r(chanindex);
      multisphere.o(i,:) = vol.o(chanindex,:);
    end
    vol = multisphere;
  end
  % if the forward model is computed using the external Neuromag toolbox,
  % we have to add a selection of the channels so that the channels
  % in the forward model correspond with those in the data.
  if isfield(vol, 'type') && strcmp(vol.type, 'neuromag')
    vol.chansel = match_str(grad.label, cfg.channel);
  end
  % if the forward model is computed using the code from Guido Nolte, we
  % have to initialize the volume model using the gradiometer coil
  % locations
  if isfield(vol, 'type') && strcmp(vol.type, 'nolte')
    % compute the surface normals for each vertex point
    if ~isfield(vol.bnd, 'nrm')
      fprintf('computing surface normals\n');
      vol.bnd.nrm = normals(vol.bnd.pnt, vol.bnd.tri);
    end
    % estimate center and radius
    [center,radius]=sphfit([vol.bnd.pnt vol.bnd.nrm]);
    % initialize the forward calculation (only if gradiometer coils are available)
    if size(grad.pnt,1)>0 
      vol.forwpar = meg_ini([vol.bnd.pnt vol.bnd.nrm], center', cfg.order, [grad.pnt grad.ori]);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the electrode positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg, 'elecfile') 
  fprintf('reading electrodes from file %s\n', cfg.elecfile);
  elec = read_sens(cfg.elecfile);
elseif isfield(cfg, 'elec') 
  fprintf('using electrodes specified in the configuration\n');
  elec = cfg.elec;
elseif isfield(data, 'elec') 
  fprintf('using electrodes specified in the data\n');
  elec = data.elec;
end

% do postprocessing of electrodes
if exist('elec', 'var')
  % update the channelselection with those present in the gradiometer array
  cfg.channel = channelselection(cfg.channel, elec.label);
  % select only electrodes present in the data
  sel        = match_str(elec.label, cfg.channel);
  fprintf('selected %d electrodes\n', length(sel)); 
  elec.pnt   = elec.pnt(sel,:);
  elec.label = elec.label(sel);
  % create a 2D projection and triangulation
  elec.prj   = elproj(elec.pnt);
  elec.tri   = delaunay(elec.prj(:,1), elec.prj(:,2));
end

% do postprocessing of volume and electrodes in case of BEM model
if exist('elec', 'var') && isfield(vol, 'bnd') && ~isfield(vol, 'tra')
  % determine boundary corresponding with skin and brain
  if ~isfield(vol, 'skin')
    vol.skin   = find_outermost_boundary(vol.bnd);
  end
  if ~isfield(vol, 'source')
    vol.source  = find_innermost_boundary(vol.bnd);
  end
  if size(vol.mat,1)~=size(vol.mat,2) && size(vol.mat,1)==length(elec.pnt)
    fprintf('electrode transfer and system matrix were already combined\n');
  else
    fprintf('projecting electrodes on skin surface\n');
    % compute linear interpolation from triangle vertices towards electrodes
    el   = project_elec(elec.pnt, vol.bnd(vol.skin).pnt, vol.bnd(vol.skin).tri);
    tra  = transfer_elec(vol.bnd(vol.skin).pnt, vol.bnd(vol.skin).tri, el);
    % construct the transfer from all vertices (also brain/skull) towards electrodes
    vol.tra = [];
    for i=1:length(vol.bnd)
      if i==vol.skin
        vol.tra = [vol.tra, tra];
      else
        vol.tra = [vol.tra, zeros(size(el,1), size(vol.bnd(i).pnt,1))];
      end
    end
    vol.tra    = sparse(vol.tra); % convert to sparse matrix to speed up multiplications
    % incorporate the transfer and the system matrix into one matrix
    % this speeds up the subsequent repeated leadfield computations
    fprintf('combining electrode transfer and system matrix\n');
    vol.mat = vol.tra * vol.mat;
    vol = rmfield(vol, 'tra');
  end
  % ensure that the model potential will be average referenced
  vol.mat = avgref(vol.mat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restructure the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('elec', 'var')
  sens = elec;
  sens.type = senstype(sens);
elseif exist('grad', 'var')
  sens = grad;
  sens.type = senstype(sens);
else
  error('cannot find electrodes or gradiometers'); 
end

% update the configuration, so that the calling function exactly knows wat was selected
cfg.channel = sens.label;


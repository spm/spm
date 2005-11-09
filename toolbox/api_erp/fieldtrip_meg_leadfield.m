function [lf] = fieldtrip_meg_leadfield(pos, grad, vol);

% MEG_LEADFIELD magnetic leadfield for a dipole in a volume conductor
% This function provides support for a single sphere model and a
% multi-sphere model, and in the future it will also compute the
% leadfield for a precomputed boundary element model. This function
% computes the leadfield for a magnetic dipole if the volume conductor
% model is empty.
% 
% [lf] = meg_leadfield(pos, grad, vol);
%
% with input arguments
%   pos		position dipole (1x3 or Nx3)
%   grad	magnetometer/gradiometer definition
%   vol	volume conductor definition

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: meg_leadfield.m,v $
% Revision 1.8  2004/09/22 10:38:24  roboos
% reordered the elseifs to avoid the Neuromag BEM model being rejected
%
% Revision 1.7  2004/09/21 11:27:57  roboos
% renamed the variable volume to vol
%
% Revision 1.6  2004/09/03 09:11:22  roboos
% incorporated Joachims change for better neuromag forward model handling
%
% Revision 1.5  2003/03/24 12:30:06  roberto
% added support for multi-sphere volume conductor model
%
% Revision 1.4  2003/03/13 13:46:19  roberto
% changed warning for magnetic dipole in vacuum
%
% Revision 1.3  2003/03/12 09:28:26  roberto
% added support for magnetic dipole with empty volume
%
% Revision 1.2  2003/03/11 14:45:36  roberto
% updated help and copyrights
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single-sphere volume conductor model
if isfield(vol, 'r') & all(size(vol.r)==[1 1])

  if isfield(grad, 'pnt1') & isfield(grad, 'pnt2')
    % this appears to be an old fashioned gradiometer definition
    % add the magnetometer positions at the upper coil
    nchan = size(grad.pnt1, 1);
    pnt = [grad.pnt1; grad.pnt2];
    ori = [grad.ori1; grad.ori2];
  else
    pnt = grad.pnt;
    ori = grad.ori;
  end

  if isfield(vol, 'o')
    % shift dipole and magnetometers to origin of sphere
    pos = pos - repmat(vol.o, size(pos,1), 1);
    pnt = pnt - repmat(vol.o, size(pnt,1), 1);
  end

  if size(pos,1)>1
    % loop over multiple dipoles
    lf = zeros(size(pnt,1),3*size(pos,1));
    for i=1:size(pos,1)
      lf(:,(3*i-2):(3*i)) = fieldtrip_meg_leadfield1(pos(i,:), pnt, ori);
    end
  else
    % only single dipole
    lf = fieldtrip_meg_leadfield1(pos, pnt, ori);
  end

  if isfield(grad, 'pnt2')
    % this appears to be a gradiometer definition
    % construct the gradient from both magnetometers
    lf = lf(1:nchan, :) - lf((nchan+1):end,:);
  end

  if isfield(grad, 'tra')
    % this appears to be the modern complex gradiometer definition
    % construct the channels from a linear combination of all magnetometers
    lf = grad.tra * lf;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multi-sphere volume conductor model
elseif isfield(vol, 'r') & ~all(size(vol.r)==[1 1])

  nspheres = length(vol.r);
  if size(vol.o, 1) ~= nspheres
    error('different number of spheres for the radius and origin')
  end

  if isfield(grad, 'pnt1') & isfield(grad, 'pnt2')
    error('oldfashoned gradiometer format not supported for multi-sphere headmodel');
  end

  if size(grad.pnt)~=nspheres
    error('number of spheres is not equal to the number of coils');
  end

  ndipoles = size(pos,1);
  lf = zeros(nspheres, 3*ndipoles);
  for chan=1:nspheres
    for dip=1:ndipoles
      % shift dipole and magnetometer coil to origin of sphere
      dippos = pos(dip,:)       - vol.o(chan,:);
      chnpos = grad.pnt(chan,:) - vol.o(chan,:);
      tmp = meg_leadfield1(dippos, chnpos, grad.ori(chan,:));
      lf(chan,(3*dip-2):(3*dip)) = tmp;
    end
  end

  if isfield(grad, 'tra')
    % this appears to be the modern complex gradiometer definition
    % construct the channels from a linear combination of all magnetometers
    lf = grad.tra * lf;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use external Neuromag toolbox for forward computation
% this requires that "megmodel" is initialized
elseif isfield(vol, 'type') & strcmp(vol.type, 'neuromag')
  ndipoles = size(pos,1);
  % compute the forward model for all channels
  tmp1 = ones(1, ndipoles);
  tmp2 = 0.01*pos';  %convert to cm
  lf = megfield([tmp2 tmp2 tmp2],[[1 0 0]'*tmp1 [0 1 0]'*tmp1 [0 0 1]'*tmp1]);
  % select only those channels from the forward model that are part of the gradiometer definition
  lf = lf(vol.chansel,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magnetic dipole instead of electric (current) dipole
% which implies that there is no volume conductor but an infinite vacuum
elseif isempty(vol)
  persistent warning_issued;
  if isempty(warning_issued)
    warning('assuming magnetic dipole in infinite vacuum');
    warning_issued = 1;
  end

  if isfield(grad, 'pnt1') & isfield(grad, 'pnt2')
    % this appears to be an old fashioned gradiometer definition
    % add the magnetometer positions at the upper coil
    nchan = size(grad.pnt1, 1);
    pnt = [grad.pnt1; grad.pnt2];
    ori = [grad.ori1; grad.ori2];
  else
    pnt = grad.pnt;
    ori = grad.ori;
  end

  if size(pos,1)>1
    % loop over multiple dipoles
    lf = zeros(size(pnt,1),3*size(pos,1));
    for i=1:size(pos,1)
      lf(:,(3*i-2):(3*i)) = magnetic_dipole(pos(i,:), pnt, ori);
    end
  else
    % only single dipole
    lf = magnetic_dipole(pos, pnt, ori);
  end

  if isfield(grad, 'pnt2')
    % this appears to be a gradiometer definition
    % construct the gradient from both magnetometers
    lf = lf(1:nchan, :) - lf((nchan+1):end,:);
  end

  if isfield(grad, 'tra')
    % this appears to be a complex gradiometer definition
    % construct the channels from a linear combination of all magnetometers
    lf = grad.tra * lf;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% realistic BEM volume conductor model
elseif isfield(vol, 'bnd')
  error('sorry, BEM volume conductor model still unsupported for MEG');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unrecognized volume conductor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  error('unrecognized volume conductor model');
end



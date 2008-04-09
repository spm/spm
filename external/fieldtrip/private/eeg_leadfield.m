function [lf] = eeg_leadfield(pos, elc, vol)

% EEG_LEADFIELD electric leadfield for a dipole in spherical or realistic
% volume conductor. This provides a wrapper around the functions EEG_LEADFIELD1
% and EEG_LEADFIELD4 for spherical models and EEG_LEADFIELDB for models based on
% the boundary element method.
%
% [lf] = eeg_leadfield(pos, elc, vol)
%
% with the input arguments
%   pos		position dipole (1x3 or Nx3)
%   elc		position electrodes
%   vol		structure defining the volume conductor model
%
% a spherical volume conductor model should have the fields
%   vol.r	radius of spheres
%   vol.c	conductivity of spheres
%   vol.o	origin of sphere (can be omitted)
%   vol.t 	constant factors for 4-sphere series expansion (can be omitted)
%
% a realistical volume conductor model should have the fields
%   vol.bnd	structure array with vertices and triangles of each boundary
%   vol.cond	conductivity of all compartments
%   vol.mat 	system matrix, which can include the electrode interpolation
%   vol.tra 	transfer matrix for interpolation from vertices to electrodes
% the BEM compartment boundaries are described by a structure array with
%   vol.bnd(i).pnt
%   vol.bnd(i).pnt
%
% See also MEG_LEADFIELD, READ_ASA_VOL

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: eeg_leadfield.m,v $
% Revision 1.5  2003/10/07 08:40:43  roberto
% moved BEM specific code into a separate function
%
% Revision 1.4  2003/07/29 11:54:02  roberto
% removed ellipsiod support (untested/unused sofar)
% renamed volume into vol
% updated help and comments
%
% Revision 1.3  2003/06/03 08:29:39  roberto
% *** empty log message ***
%
% Revision 1.2  2003/03/11 14:45:36  roberto
% updated help and copyrights
%

% approximate timing per EEG channel on a PIII/800Hz
%              m-code   p-code  mex-file
% 1 sphere      0.30     0.30     0.15
% 4 sphere      0.56     0.56     1.84 (!!)            with Nmax=10
% 4 sphere      1.46     1.46     5.36 (!!)            with Nmax=30

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spherical volume conductor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(vol, 'r')
  
  % FIXME, this is not consistent between spherical and BEM
  % sort the spheres from the smallest to the largest
  [vol.r, indx] = sort(vol.r);
  vol.c = vol.c(indx);
  
  Nspheres = length(vol.c);
  if length(vol.r)~=Nspheres
    error('the number of spheres in the volume conductor model is ambiguous');
  end
  
  if isfield(vol, 'o')
    % shift the origin of the spheres, electrodes and dipole
    elc = elc - repmat(vol.o, size(elc,1), 1);
    pos = pos - repmat(vol.o, size(pos,1), 1);
  end

  if Nspheres==1
    if size(pos,1)>1
      % loop over multiple dipoles
      lf = zeros(size(elc,1),3*size(pos,1));
      for i=1:size(pos,1)
        lf(:,(3*i-2):(3*i)) = eeg_leadfield1(pos(i,:), elc, vol);
      end
    else
      % only single dipole
      lf = eeg_leadfield1(pos, elc, vol);
    end

  elseif Nspheres==2
    vol.r = [vol.r(1) vol.r(2) vol.r(2) vol.r(2)];
    vol.c = [vol.c(1) vol.c(2) vol.c(2) vol.c(2)];
    if size(pos,1)>1
      % loop over multiple dipoles
      lf = zeros(size(elc,1),3*size(pos,1));
      for i=1:size(pos,1)
        lf(:,(3*i-2):(3*i)) = eeg_leadfield4(pos(i,:), elc, vol);
      end
    else
      % only single dipole
      lf = eeg_leadfield4(pos, elc, vol);
    end

  elseif Nspheres==3
    vol.r = [vol.r(1) vol.r(2) vol.r(3) vol.r(3)];
    vol.c = [vol.c(1) vol.c(2) vol.c(3) vol.c(3)];
    if size(pos,1)>1
      % loop over multiple dipoles
      lf = zeros(size(elc,1),3*size(pos,1));
      for i=1:size(pos,1)
        lf(:,(3*i-2):(3*i)) = eeg_leadfield4(pos(i,:), elc, vol);
      end
    else
      % only single dipole
      lf = eeg_leadfield4(pos, elc, vol);
    end

  elseif Nspheres==4
    if size(pos,1)>1
      % loop over multiple dipoles
      lf = zeros(size(elc,1),3*size(pos,1));
      for i=1:size(pos,1)
        lf(:,(3*i-2):(3*i)) = eeg_leadfield4(pos(i,:), elc, vol);
      end
    else
      % only single dipole
      lf = eeg_leadfield4(pos, elc, vol);
    end

  else
    error('more than 4 concentric spheres are not supported');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% realistic BEM volume conductor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isfield(vol, 'bnd')

  lf = eeg_leadfieldb(pos, elc, vol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unrecognized volume conductor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  error('unrecognized volume conductor model');
end


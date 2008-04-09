function [lf] = eeg_leadfieldb(pos, elc, vol)

% EEG_LEADFIELDB computes the electric leadfield for a dipole in a volume
% using the boundary element method
%
% [lf] = eeg_leadfieldb(pos, elc, vol)
%
% with the input arguments
%   pos		position dipole (1x3 or Nx3)
%   elc		position electrodes (optional, can be empty)
%   vol		volume conductor model
%
% the volume conductor model is a structure and should have the fields
%   vol.bnd	structure array with vertices and triangles of each boundary
%   vol.cond	conductivity of all compartments
%   vol.mat 	system matrix, which can include the electrode interpolation
%   vol.tra 	transfer matrix for interpolation from vertices to electrodes
%
% the compartment boundaries are described by a structure array with
%   vol.bnd(i).pnt
%   vol.bnd(i).pnt

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: eeg_leadfieldb.m,v $
% Revision 1.2  2005/12/06 11:41:07  roboos
% added support for dipoli models
% restructured the whole code
% combine electrode transfer and system matrix in a single matrix to speed forward computations up
%
% Revision 1.1  2003/10/07 08:40:22  roberto
% made separate function for BEM, based on part of eeg_leadfield.m
%

isavo    = 0;  % system matrix computed using code from Adriaan van Oosterom
isdipoli = 0;  % system matrix computed using Thom Oostendorp's DIPOLI
isasa    = 0;  % system matrix computed using ASA from www.ant-neuro.com

% determine the type of BEM system matrix
if isfield(vol, 'type') && strcmp(vol.type, 'avo')
  isavo = 1;
elseif isfield(vol, 'type') && strcmp(vol.type, 'dipoli')
  isdipoli = 1;
else
  isasa = 1;
end

% determine the number of compartments
ncmp = length(vol.cond);

% do some tests for the correctness of the model
if ~isfield(vol, 'bnd')
  error('there are no compartment boundaries present');
end

if length(vol.bnd)~=ncmp
  error('the number of compartments in the volume in ambiguous');
end

if ~isfield(vol, 'mat')
  error('there is no BEM system matrix present');
end

if isfield(vol, 'tra') & size(vol.tra,2)~=size(vol.mat,1)
  error('ambiguity between electrode transfer matrix and system matrix');
end

if isfield(vol, 'tra') & ~isempty(elc) & size(vol.tra,1)~=size(elc,1)
  error('ambiguity between number of electrodes and transfer matrix');
end

if ~isfield(vol, 'tra') & ~isempty(elc) & size(vol.mat,1)~=size(elc,1)
  error('ambiguity between number of electrodes and system matrix');
end

if ~isfield(vol, 'source')
  % determine the source compartment
  if ncmp==1
    % there is only one compartment
    vol.source = 1;
  else
    % try to locate the innermost compartment and
    % assume that the sources are in the innermost compartment
    innermost = find_innermost_boundary(vol.bnd);
    warning(sprintf('assuming that BEM compartment %d contains the sources', innermost));
    vol.source = innermost;
  end
end
% determine the conductivity of the source compartment
cond = vol.cond(vol.source);

% compute the infinite medium potential on all vertices
if isavo
  % the code by Adriaan van Oosterom does not implement isolated source approach
  lf = [];
  for i=1:ncmp
    lf = [lf; inf_medium_leadfield(pos, vol.bnd(i).pnt, mean(vol.sigmas(i,:)))];
  end
elseif isdipoli
  % concatenate the vertices of all compartment boundaries in a single Nx3 matrix
  pnt = [];
  for i=1:ncmp
    pnt = [pnt; vol.bnd(i).pnt];
  end
  % dipoli incorporates the conductivity into the system matrix
  lf = inf_medium_leadfield(pos, pnt, 1);
elseif isasa
  % concatenate the vertices of all compartment boundaries in a single Nx3 matrix
  pnt = [];
  for i=1:ncmp
    pnt = [pnt; vol.bnd(i).pnt];
  end
  % assume that isolated potential approach was used
  lf = inf_medium_leadfield(pos, pnt, cond);
end

% compute bounded medium potential on all vertices
lf = vol.mat * lf;

% optionally do bilinear interpolation from vertices towards electrodes
if isfield(vol, 'tra')
  lf = vol.tra * lf;
end

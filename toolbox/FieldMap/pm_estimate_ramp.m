function varargout = pm_estimate_ramp(varargin)
% Estimate and remove linear phase-ramps in the x-, y- and z-direction in pm
% FORMAT [ramps,pm] = pm_estimate_ramp(pm,mask)
%
% Input:
% pm     : 2 or 3D phasemap that has not been unwrapped
% mask   : Mask that indicates which voxels are worth
%          bothering with and which are not.
%
% Output:
% ramps  : 3x1 vector with ramp values (radians/voxel) in the x-,
%          y- and z-direction. Estimated as mean(wrap(phi(x+1)-phi(x)))
%          for the x-direction and equivalently for the other directions.
% pm     : Same as pm in, but with linear ramps removed.
%
% This routine was written on the suggestion of Mark J, and will
% potentially improve performance of subsequent phase-unwrapping. I haven't
% actually found it particularly helpful, and it may simply have been a
% sneaky fMRIB attempt to delay the SPM phasemap toolbox.
%__________________________________________________________________________

% Jesper Andersson
% Copyright (C) 2003-2022 Wellcome Centre for Human Neuroimaging

error('mex-function pm_estimate_ramp.c not compiled');

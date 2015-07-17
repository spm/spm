function varargout = pm_pad(varargin)
% Pads a (partially) unwrapped phasemap such that the phase
% at a non-unwrapped location is a weighted average of unwrapped
% neighbouring phase-values.
% FORMAT [pm,wmap] = pm_pad(pm,wmap,kernel)
%
% Input:
% pm     : 2 or 3D phasemap where some voxels have been unwrapped 
%          and some not.
% wmap   : Wrap-map, where a non-zero value indicates corresponding 
%          phase-value in pm has been unwrapped.
% kernel : kernel used to generate a weighted average of surrounding
%          voxels.
%
% Output: 
% pm     : Same as pm in, but where some previously unwrapped
%          phase-values have now been replaced.
% wmap   : Same as wmap in, but where values that was replaced
%          by weighted average in pm have now been set.
%__________________________________________________________________________
% Copyright (C) 2008-2015 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_pad.m 6501 2015-07-17 14:32:09Z spm $

error('mex-function pm_pad.c not compiled');

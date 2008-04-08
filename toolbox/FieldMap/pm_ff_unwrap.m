function varargout = pm_ff_unwrap(varargin)
% Performs phase-unwrapping along the shore of a flood-fill
% progressing from some seed in a low-variance area towards
% higher variance areas.
% FORMAT [pm,wm] = pm_ff_unwrap(pm,vm,wm,mask,thres)
%
% Input:
% pm     : 2 or 3D phasemap where some voxels have been unwrapped 
%          and some not.
% vm     : Variance map, indicating the variance of the phase estimate.
% wm     : Wrap-map, where a non-zero value indicates corresponding 
%          phase-value in pm has been unwrapped.
% mask   : Mask that indicates which voxels are worth
%          bothering with and which are not.
% thres  : Array of thresholds such that unwrapping proceeds by
%          serial dilates from some seed point in wm until all
%          connected voxels with vm<thres(1) have been unwrapped,
%          after which the same happens with vm<thres(2) etc etc
%
% Output:
% pm     : Same as pm but with unwrapping now completed as far as
%          indicated by wm on output.
% wm     : Indicates how far unwrapping has progressed.
%
% Note that the following are equivalent in terms of the final
% result
%
% [pm,wm] = pm_ff_unwrap(pm,vm,wm,mask,thres);
%
% for i=1:length(thres)
%   [pm,vm] = pm_ff_unwrap(pm,vm,wm,maskm,thres(i));
% end
%
% where the later would be slightly slower but would allow
% reporting of progress to some GUI.
%_________________________________________________________________
% Jesper Andersson 30/9-03
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_ff_unwrap.m 1317 2008-04-08 16:16:38Z chloe $

error('mex-function pm_ff_unwrap.c not compiled');



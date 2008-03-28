function varargout = spm_bias_mex(varargin)
% A mex function involved in bias correction
%
% FORMAT [Alpha,Beta,ll, h, n] = spm_bias_mex(V,B1,B2,B3,T,[mx nh])
%   V           - a handle (spm_vol) to an image.
%   B1, B2 & B3 - basis functions in x, y and z directions
%   T           - basis function coefficients
%   mx          - maximum intensity in image
%   nh          - number of bins in histogram
%
%   Alpha       - Matrix of second derivatives of ll w.r.t. T
%   Beta        - Vector of first derivatives of ll w.r.t. T
%   ll          - log-likelihood
%   h           - Histogram of bias corrected data
%   n           - Number of voxels
%
% This mex function is called by spm_bias_estimate.m as part of a non-parametric
% bias correction algorithm. The mex function computes a histogram of the bias
% corrected (according to latest estimates of coefficients) image in order to
% compute the log-likelihood and its first and second derivatives.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging


% John Ashburner
% $Id: spm_bias_mex.m 1271 2008-03-28 15:06:48Z john $

[pth,nam,ext ] = fileparts(mfilename);
error('The function "%s" is not compiled for %s in MATLAB %s.\nSee %s%csrc%cMakefile for information about how you may be able to compile it.\n',...
      nam, computer,version,spm('dir'),filesep,filesep);


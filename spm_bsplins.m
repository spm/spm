function varargout = spm_bsplins(varargin)
% Sample a volume using B-spline interpolation
% FORMAT [f,dfx,dfy,dfz] = spm_bsplins(c,x,y,z,o)
% 	c - volume of B-spline coefficients (from spm_bsplinc)
% 	x,y,z - co-ordinates of sampled points
% 	o - order of B-splines (must be same as used by spm_bsplinc)
% 	f - sampled data
% 	dfx,dfy,dfz - sampled first derivatives
%
% This function takes B-spline basis coefficients from spm_bsplinc,
% and re-convolves them with B-splines centred at the new sample points.
%
%_______________________________________________________________________
%
% References:
%	M. Unser, A. Aldroubi and M. Eden.
%	"B-Spline Signal Processing: Part I-Theory,"
%	IEEE Transactions on Signal Processing 41(2):821-832 (1993).
%
%	M. Unser, A. Aldroubi and M. Eden.
%	"B-Spline Signal Processing: Part II-Efficient Design and
%	Applications,"
%	IEEE Transactions on Signal Processing 41(2):834-848 (1993).
%
%	M. Unser.
%	"Splines: A Perfect Fit for Signal and Image Processing,"
%	IEEE Signal Processing Magazine, 16(6):22-38 (1999)
%
%	P. Thévenaz and T. Blu and M. Unser.
%	"Interpolation Revisited"
%	IEEE Transactions on Medical Imaging 19(7):739-758 (2000).
%_______________________________________________________________________
% %W% John Ashburner %E%

%-This is merely the help file for the compiled routine
error('spm_bsplins.c not compiled.');


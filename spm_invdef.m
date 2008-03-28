function varargout = spm_invdef(varargin)
% Create the inverse of a deformation field.
% FORMAT [Y1,Y2,Y3] = spm_invdef(X1,X2,X3,dimY,MX,MY)
%   X1, X2 & X3 - Three components of original deformation field.
%                 Note that these point from voxels to a coordinate
%                 system in mm.
%   dimY        - A three element vector encoding the dimensions of
%                 the inverse of the deformation.
%   MX          - The voxel-to-world mapping for the elements of the
%                 original deformation.
%   MY          - The voxel-to-world mapping for the elements of the
%                 resulting inverse deformation.
%
%   Y1, Y2 & Y3 - Three components of inverse deformation field.
%                 Note that these point from voxels to a coordinate
%                 system in mm.
%
% Deformations are encoded as piecewise affine transforms. The space
% between each set of 8 neighbouring voxels is divided into 5
% tetrahedra, where there is an affine mapping within each of them.
%_______________________________________________________________________
% Inverting the deformation is as described in the appendix of:
%
%     John Ashburner, Jesper L.R. Andersson and Karl J. Friston.
%     "Image Registration Using a Symmetric Prior in Three Dimensions".
%     Human Brain Mapping 9:212-225(2000)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_invdef.m 1271 2008-03-28 15:06:48Z john $

[pth,nam,ext ] = fileparts(mfilename);
error('The function "%s" is not compiled for %s in MATLAB %s.\nSee %s%csrc%cMakefile for information about how you may be able to compile it.\n',...
      nam, computer,version,spm('dir'),filesep,filesep);


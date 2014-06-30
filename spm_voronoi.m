function [vor, dist] = spm_voronoi(img, seeds, distance)
% Geodesic Discrete Voronoi Diagram - a compiled routine
% FORMAT [vor, dist] = spm_voronoi(img, seeds, distance)
%
% img       - binary image:  > 0  : inside
%                            <= 0 : outside 
% seeds     - {n x 3} array of the n seeds positions [in voxels]
% distance  - type of chamfer distance to use ('d4', 'd8', 'd34' or 'd5711')
%             (default is 'd34')
%
% vor       - Geodesic Discrete Voronoi diagram 
%             (label is equal to the index of the seed in 'seeds')
% dist      - Geodesic Distance map of img with seeds as objects
%
% Compute the geodesic discrete Voronoi Diagram of an image of labelled 
% objects using front propagation. The distance map is also available 
% on output.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_voronoi.m 6079 2014-06-30 18:25:37Z spm $

%-This is merely the help file for the compiled routine
error('spm_voronoi.c not compiled - see Makefile')

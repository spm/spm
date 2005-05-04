% Create 2D scatter-plot of two images after affine transformation
% FORMAT H = spm_hist2(G,F,M,s)
% G  - unsigned 8 bit 3D array representing the first volume
% F  - unsigned 8 bit 3D array representing the second volume
% M  - the affine transformation matrix so that G(x) is plotted
%      against F(x)
% s  - 3 element vector representing the sampling density (in voxels).
%      A value of [1 1 1] means that approximately all voxels are
%      sampled, whereas [4 4 4] means that only 1 in 64 voxels are
%      sampled.
%
% This function is called by spm_mireg for rapidly computing
% joint histograms for mutual information image registration.
%
% Note that the function jitters the sampling of the data to reduce
% interpolation artifacts.
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_hist2.m 112 2005-05-04 18:20:52Z john $


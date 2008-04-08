function pm = pm_merge_regions(pm,rima,ii,jj,nn,pp,rs)
% Merges regions as defined in connectogram to minimise 
% total costfunction (sum of phase-differences across 
% region borders).
% FORMAT: pm = pm_merge_regions(pm,rima,ii,jj,nn,pp,rs)
%
% Input:
% pm       : Phase-map
% rima     : Label map consisting of connected regions indentified
%            by unique labels. Use pm_initial_regions to get rima.
% ii       : Array of row indicies.
% jj       : Array of column indicies.
% nn       : Array of no. of voxels in borders between regions.
%            So e.g. if ii[10]=5, jj[10]=9 and nn[10]=123 it 
%            means that regions 5 and 9 have a common border
%            (are connected) and that this border has 123 voxels. 
% pp       : Array of sum of phase differences between regions.
%            So e.g. if ii[10]=5, jj[10]=9 and pp[10]=770.2 it 
%            means that regions 5 and 9 have a common border
%            (are connected) and that for paired voxels across
%            this border the sum of phase differenes is 770.2.
%            N.B. the subtraction is phi(ii(i))-phi(jj(i)), 
%            which in the example above means that the phase is
%            smaller in region 9 than in region 5.
% rs       : List of region sizes, so that e.g. if rs[13]=143 it
%            means that the regions with label 13 consists 
%            of 143 voxels.
%
% Output:
% pm       : Phase-map after merging of all regions in rima that
%            are connected.
%
% This routine is based on the MRM paper by Mark J. Very briefly it will
% use the summary statistic in the matrices N and P, where each entry in
% N signifies the number of voxels along the common border of the regions
% whose labels correspond to row and column of the matrix. E.g. N(i,j) (for i<j) 
% signifies the number of voxels along the border between regions labelled
% i and j. The matrix P is organised in the same manner, with the difference
% that the numbers correspond to the sum of differences of phase values
% across voxel-faces along that border. The direction of the difference has
% been (arbitrarily) chosen such that we take phi(i)-phi(j) where i<j.
%
% Now we want to merge all these regions, such that after merging all
% phase-wraps will have been resolved. An assumption here is that any
% phase-wraps will always be along borders of the initial regions,
% something that is (almost) guaranteed by the way in which we create them.
% 
% There are two aspects to the merging
% 1. We want to detect and correct for any phase-wrap between regions
%    i and j when merging them.
% 2. We want to merge the regions in such an order that more "important"
%    regions are merged first. This is functionally similar to the 
%    progression of wrapping from low->high varinace areas in region-growing
%    approches.
%
% The first goal is easily reached by noting that (P(i,j)/N(i,j))/2pi
% is a good guess for the number of wraps that differ between regions
% i and j.
%
% The second goal is reached by merging the pairs of regions that have
% the largest border (i.e. the largest N(i,j)) first (it is a little
% more elaborate, but basically like that). 
%
% The rest is really just about being really careful when updating the
% stats regarding all the connections between a newly merged regions
% and all the regions that bordered to one or both of the regions
% constituting the new region.
%
% Jenkinson M. 2003. Fast, automated, N-dimensional phase-unwrapping 
% algorithm. MRM 49:193-197.
%_________________________________________________________________________
% Jesper Andersson 2/10-03
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_merge_regions.m 1317 2008-04-08 16:16:38Z chloe $

error('pm_merge_regions.c has not been compiled');

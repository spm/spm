function seed = pm_seed(angvar,mask,pxs)
% Find a suitable (hopefully) seed point from which
% to start watershed-based unwrapping.
% FORMAT: seed = pm_seed(angvar,mask,pxs)
%
% Input:
% angvar  : Map of variance of (voxelwise) estimates
%           of phase angle.
% mask    : Tells us which part of angvar to consider.
% pxs     : Array of voxel sizes, used to ensure
%           isotropic smoothing.
%
% Output:
% seed    : Coordinates of suitable seed point.
%
% In order to find a seed point we first threshold the
% variance map at a quarter of the variance of a U(-pi,pi)
% distribution. This gives us a binary image with ones only
% for low variance regions. This is then smoothed with a
% very wide gaussian kernel (50mm). The maximum of
% the smoothed map is then pretty much a centre-of-mass
% of the "low-variance volume". It could however in 
% principle be a relatively high variance voxel 
% surrounded by low-variance voxels. Therefore we pick
% a percentage of the highest voxels in the smooth map
% (i.e. we pick a neighbourhood) and then pick the location
% of those that has the lowest variance in the original
% variance map.
%___________________________________________________________
% Jesper Andersson 1/10-03 
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_seed.m 4842 2012-08-15 18:02:30Z guillaume $

if nargin < 3
   mask = ones(size(angvar));
end

%
% First let us create a volume where a high
% value indicates a "high density of relatively
% low variance voxels".
%
dim = size(angvar);
if length(dim) == 2 dim(3) = 1; end

M = eye(4)*diag([pxs(1) pxs(2) pxs(3) 1]);
M = M - [zeros(4,3) M*[mean(1:dim(1)) mean(1:dim(2)) mean(1:dim(3)) 0]']; 

P = struct('dim',     [dim 64],...
           'pinfo',   [1 0]',...
           'mat',     M);
P.dat = double(angvar<(pi^2)/12);
svol = zeros(size(angvar));
spm_smooth(P.dat,svol,50);

%
% A high value in svol "probably" indicates a
% voxel with a low variance, surrounded by a 
% lot of other voxels with low variance (i.e.
% a good place to start unwrapping from). 
% However, it COULD also be a voxel with
% "not so low variance" surrounded by low
% variance voxels. To avoid that trap we pick
% the voxel with the lowest variance in the
% unsmoothed variance map, out of the 5% of
% voxels in the svol with the highest values.
% It's all very heuristic and ugly.
%

[N,X] = hist(svol(logical(mask(:))),100);
indx = find(cumsum(N)>0.95*length(find(mask)));
thres = X(indx(1)-1);

indx = find(svol(:)>thres);
[mv,mi] = min(angvar(indx));
seed = zeros(1,3);
[seed(1),seed(2),seed(3)] = ind2sub(size(angvar),indx(mi));

return

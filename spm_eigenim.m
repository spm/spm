function [e,l,E] = spm_eigenim(P,mc,voxels,n,sigma,GM)
% Eigenimages and eigenvectors from image files
% FORMAT [e,l,E] = spm_eigenim(P,mc,voxels,n,sigma,GM)
%
% P        - list of file-names or -headers
% mc       - Mean correct, 0=not, 1=rowwise, 2=columnwise, 3=both
% voxels   - List of voxels (xyz) to use
% n        - Number of eigenimages to compute
% sigma    - Size of filter kernel. 6/sqrt(8*log(2))/RT is sensible
%__________________________________________________________________________
%
% This routine is used to calculate the n first eigenimages for the
% image series pointed to by P. In order to conserve memory the SVD is
% not used, but rather the eigenvalue decomposition of the scan-by-scan
% covariance matrix. The covariance matrix is calculated in "steps" by
% reading an appropriate amount of voxels at a time, calculating the
% partial covariance matrices and adding them up.
%__________________________________________________________________________
% %W% Jesper Andersson %E%


if (~size(P)) % If no files given
   P = spm_get(+Inf,'.img','Select data for eigenimage analysis');
end

% Get header data if not already there

if ~isfield(P,'dim') P = spm_vol(P); end

% Get point list if no voxels given.

if ~size(voxels)
   voxels = zeros(prod(P(1).dim(1:3)),3);
   [x,y,z] = meshgrid(1:P(1).dim(1),1:P(1).dim(2),1:P(1).dim(3));
   voxels(:,1) = x(:);
   voxels(:,2) = y(:);
   voxels(:,3) = z(:);
end
clear x y z;

% If correction for global mean should be done
% and globals not being passed to routine.

if (cm==2 | cm==3) & ~size(GM)
   GM = zeros(1,size(P,1));
   for i=1:size(P,1)
      vol = spm_sample_vol(P(i),voxels(:,1),voxels(:,2),voxels(:,3),0);
      GM(i) = mean(vol);
   end
end
   

if ~size(sigma) sigma = 0; 
else K = spm_sptop(sigma,size(P,1)); end
   
% Determine chunk size aiming at ~8MB per chunk

cz = round(1e6 / size(P,1));

% Create scan-by-scan covariance matrix

cm = zeros(size(P,1));
X = zeros(cz,size(P,1));
for i=1:cz:size(voxels,1)
   if i+cz > size(voxels,1) cz = size(voxels,1)-i; end 
   for j=1:size(P,1)
      X(1:cz,j) = spm_sample_vol(P(j),voxels(i:i+cz-1,1),voxels(i:i+cz-1,2),voxels(i:i+cz-1,3),0);
   end
   if sigma X = X*K'; end
   if mc==1 | mc==3 X(1:cz,:) = X(1:cz,:) - repmat(mean(X(1:cz,:)')',1,size(P,1)); end
   if mc==2 | mc==3 X(1:cz,:) = X(1:cz,:) - repmat(GM,cz,1); end
   cm = cm + spm_atranspa(X(1:cz,:));
end
clear X;

% Get eigenvectors of scan-by-scan matrix

[e,l] = eig(cm);

% Sort it in descending order of eigenvalues

[sl,indx] = sort(diag(l));
l = flipud(fliplr(l([indx],[indx])));
e = fliplr(e(:,[indx]));

% Because of the mean correction the matrix
% is rank deficient, why the last eigenvalue
% should be zero, and may because of round off
% errors be negative.

e(size(e,1),size(e,2)) = 0;

% Get n eigenvectors of voxel-by-voxel matrix
% if asked for it.

if nargin > 2
   if ~size(n) n = size(P,1); end
   el = e*inv(sqrt(l));
   el = el(:,1:n);
   E = zeros(size(voxels,1),n);
   xvec = zeros(size(voxels,1),1);
   for i=1:size(P)
      xvec = spm_sample_vol(P(i),voxels(:,1),voxels(:,2),voxels(:,3),0);
      E = E + xvec*el(i,:);
   end
end

return;

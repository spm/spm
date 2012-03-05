function spm_orientations(P)
% Show the orientations that SPM assumes that the data are
% stored in.  Standard Analyze format axial images will
% normally be reported as 'RPI   Left-handed'.  Some people
% will represent their axial images as Right-handed.
% 'RPI' means that the fastest changing direction (i.e.
% the first element of the voxel coordinate) in the
% file is Right->left, the middle (second element of
% voxel coordinate) is Posterior->anterior and the
% slowest (third element - indicating slice number) is
% Inferior->superior.
%
% One thing to watch out for is the image orientation. The
% proper Analyze format uses a left-handed co-ordinate system,
% whereas Talairach uses a right-handed one. In SPM99, images
% were flipped at the spatial normalisation stage (from one
% co-ordinate system to the other). In SPM2, a different
% approach is used, so that either a left- or right-handed
% co-ordinate system is used throughout. The SPM2 program is
% told about the handedness that the images are stored with by
% the spm_flip_analyze_images.m function and the
% defaults.analyze.flip parameter that is specified in the
% spm_defaults.m file. These files are intended to be
% customised for each site. If you previously used SPM99 and
% your images were flipped during spatial normalisation, then
% set defaults.analyze.flip=1. If no flipping took place, then
% set defaults.analyze.flip=0.
%
% Check that when using the Display facility (possibly after
% specifying some rigid-body rotations) that:
%
%     * The top-left image is coronal with the top (superior)
%       of the head displayed at the top and the left shown on
%       the left. This is as if the subject is viewed from
%       behind.
%
%     * The bottom-left image is axial with the front
%       (anterior) of the head at the top and the left shown
%       on the left. This is as if the subject is viewed from
%       above.
%
%     * The top-right image is sagittal with the front
%       (anterior) of the head at the left and the top of the
%       head shown at the top. This is as if the subject is
%       viewed from the left.

%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_orientations.m 4678 2012-03-05 18:01:33Z john $

if nargin<1,
    P = spm_select(Inf,'image','Select the images...');
end

if spm_flip_analyze_images
    fprintf('SPM is assuming left-handed storage when handedness is not indicated by the .hdr or .mat (flip=1)\n');
else
    fprintf('SPM is assuming right-handed storage when handedness is not indicated by the .hdr or .mat (flip=0)\n');
end

for i=1:size(P,1),
    Nii = nifti(P(i,:));
    M   = Nii.mat;
    [U,S,V] = svd(M(1:3,1:3));
    M   = U*V';
    lab = 'LRPAIS';
    d   = [1 -1  0  0  0  0
         0  0  1 -1  0  0
         0  0  0  0  1 -1];
    dp = M\d;
    c  = '   ';
    for j=1:3,
        [unused,ind] = max(dp(j,:));
        c(j)         = lab(ind);
    end;

    if det(M)>0,
        h = 'Right';
    else
        h = ' Left';
    end;
    fprintf('%s  %s-handed  %s\n',c,h,P(i,:));
end

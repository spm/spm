function spm_depth(P)
% Create an image of cortical depth.
% FORMAT spm_depth(P)
% P - filenames for GM and WM maps.
%
% A skeletonisation is saved, along with a map of
% cortical depth.
%
% In principle, cortical depth can be determined from maps
% of GM and WM. However, opposite banks of sulci usually touch
% making it difficult to identify the outer surface of the brain.
%
% This code uses maps of GM and WM to:
% 1. Obtain a map of background (i.e. BG = 1-GM-WM).
% 2. Takes a binary map of GM+BG>0.5, and does a topology
%    preserving conditional erosion (conditional on BG>0.5)
%    to try to identify voxels running along the medial surface
%    of sulci. The end result is a sort of binary skeletonisation
%    that tries to generate medial surfaces.
% 3. Combine the GM, WM and BG with the skeletonisation, giving
%    mGM, mWM and mBG.
% 4. Compute Euclidean distances of all voxels from the mWM map (d2).
% 5. Compute Euclidean distances of all voxels from the mBG map (d3).
% 6. Save a depth map computed from d3./(d2+d3)
%
% The end result is not perfect because it is susceptible to segmentation
% imperfections, particularly in relation to isolated speckles of WM, or
% incorrect WM topology. Also, the medial surface is mid-way between the
% inner and outer surfaces, which may not always be optimal. Thresholding
% segmentation probabilities at 0.5 is a bit arbitrary, but a probabilistic
% topology-preserving erosion would be very slow.
% 
% While the results may not be exactly correct, they are probably a
% better approximation than would be obtained from not doing the
% skeletonisation. Cortical regions to worry about most will be those
% close to the cerebellum.  It will also do strange things for cortical
% regions close to subcortical GM structures.
%
% 

if nargin<1
    P = spm_select(2,'^c.*\.nii','Select a GM & WM map');
end

% Read the GM and WM, and compute BG
Nii = nifti(P);
p1  = squeeze(Nii(1).dat(:,:,:));
p2  = squeeze(Nii(2).dat(:,:,:));
p3  = max(1-p1-p2,0);


fprintf('\nThinning\n');
skel = spm_thin(p1+p3>0.5,p3>0.5);
%skel = convn(skel,ones(3,3,3)/27,'same');

% Write out the skeleton for debugging purposes
Nio            = Nii(1);
[dr,nam,ext]   = fileparts(Nio.dat.fname);
Nio.dat.fname  = fullfile(dr,['skel_' nam '.nii']);
Nio.descrip    = 'Skeletonisation';
create(Nio);
Nio.dat(:,:,:) = skel;

% Combine the skeletonisation with the GM, WM and BG.
p3  = max(p3,skel);
p2  = min(p2,1-skel);
p1  = min(p1,1-skel);

% Voxel sizes
vx = sqrt(sum(Nii(1).mat(1:3,1:3).^2));

fprintf('Dist from WM: ');
d2 = spm_distance3(p2,vx,10);

fprintf('Dist from CSF:');
d3 = spm_distance3(p3,vx,10);

% Depth from outer surface (0) to inner surface (1)
d  = d3./(d2+d3);

% Write the output
Nio            = Nii(1);
[dr,nam,ext]   = fileparts(Nio.dat.fname);
Nio.dat.fname  = fullfile(dr,['depth_' nam '.nii']);
Nio.descrip    = 'Depth';
create(Nio);
Nio.dat(:,:,:) = d;


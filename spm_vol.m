function V = spm_vol(P)
% Get header information etc for images.
% FORMAT V = spm_vol(P)
% P - a matrix of filenames.
% V - a vector of structures containing image volume information.
% The elements of the structures are:
%       V.fname - the filename of the image.
%       V.dim   - the x, y and z dimensions of the volume, and the
%                 datatype of the image.
%       V.mat   - a 4x4 affine transformation matrix mapping from
%                 voxel coordinates to real world coordinates.
%       V.pinfo - plane info for each plane of the volume.
%              V.pinfo(1,:) - scale for each plane
%              V.pinfo(2,:) - offset for each plane
%                 The true voxel intensities of the jth image are given
%                 by: val*V.pinfo(1,j) + V.pinfo(2,j)
%              V.pinfo(3,:) - offset into image (in bytes).
%                 If the size of pinfo is 3x1, then the volume is assumed
%                 to be contiguous and each plane has the same scalefactor
%                 and offset.
%____________________________________________________________________________
%
% The fields listed above are essential for the mex routines, but other
% fields can also be incorporated into the structure.
% Unlike with spm_map_vol, the images are not memory mapped at this step,
% but are mapped when the mex routines using the volume information are
% called.
% This is a replacement for the spm_map_vol and spm_unmap_vol stuff
% that I expect will become redundant.
%_______________________________________________________________________
% %W% John Ashburner %E%

for i=1:size(P,1)
	p        = deblank(P(i,:));

	[dim vox scale dtype offset] = spm_hread(p);
	mat = spm_get_space(p);

	V(i) = struct('fname',p,'dim',[dim dtype],'mat',mat,'pinfo',[scale 0 offset]');
end

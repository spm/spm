% Get the best guess for the space of the image.
% FORMAT M = spm_get_space(imagename)
%
% imagename - the name of the image
% M         - a transformation matrix which maps voxel coordinates
%             to mm with referance to some origin.
%             The corner voxels of the image are at 1,1,1.
%___________________________________________________________________________
%
% If imagename.mat exists, and this contains matrix "M", then this is used.
% Otherwise, the image is assumed to be transverse, with the origin of the
% space found in the ORIGIN field of the header.
% If this ORIGIN field is set to [0 0 0], then the origin is assumed to be
% at the centre (center) of the volume.
%

% %W% John Ashburner FIL %E%

function M = spm_get_space(imagename)

% If we can read M from a ".mat" file with the corresponding
% name, then use this information.
matname = [spm_str_manip(imagename,'sd') '.mat'];
fp = fopen(matname,'r');
if (fp ~= -1)
	fclose(fp);
	load(matname);
	if (exist('M') == 1)
		return;
	end
end

% Assume transverse slices.

[dim vox scale typ offset origin descrip] = spm_hread([spm_str_manip(imagename,'sd') '.img']);

% If origin hasn't been set, then assume
% it is the centre of the image.
if (all(origin == 0))
	origin = dim(1:3)/2;
end

if (all(vox == 0))
	vox = [1 1 1];
end

off = -vox.*origin;
M = [vox(1) 0 0 off(1) ; 0 vox(2) 0 off(2) ; 0 0 vox(3) off(3) ; 0 0 0 1];

function M = spm_get_space(imagename, mat)
% Get or set the best guess for the space of the image.
% FORMAT M = spm_get_space(imagename) OR spm_get_space(imagename,mat)
%
% imagename - the name of the image
% M         - a transformation matrix which maps voxel coordinates
%             to mm with referance to some origin.
%             The corner voxels of the image are at 1,1,1.
% mat       - a transformation matrix.
%___________________________________________________________________________
%
% If imagename.mat exists, and this contains matrix "M", then this is used.
% Otherwise, the image is assumed to be transverse, with the origin of the
% space found in the ORIGIN field of the header.
% If this ORIGIN field is set to [0 0 0], then the origin is assumed to be
% at the centre (center) of the volume.
%___________________________________________________________________________
% %W% John Ashburner FIL %E%

imagename = deblank(imagename);
matname = [spm_str_manip(imagename,'sd') '.mat'];
if (nargin == 1)
	% If we can read M from a ".mat" file with the corresponding
	% name, then use this information.
	if (exist(matname) == 2)
		load(matname);
		if (exist('M') == 1)
			return;
		end
	end

	if exist([spm_str_manip(imagename,'sd') '.hdr']) == 2,
		% Read as much information from the ANALYZE header.
		% Assume transverse slices.
		[dim vox scale typ offset origin descrip] = spm_hread([spm_str_manip(imagename,'sd') '.img']);
		if isempty(dim),
			error(sprintf('Can''t read header for "%s"\n', [spm_str_manip(imagename,'sd') '.img']));
		end

		% If origin hasn't been set, then assume
		% it is the centre of the image.
		if all(origin == 0), origin = (dim(1:3)+1)/2; end;
		if all(vox == 0), vox = [1 1 1]; end;
		off = -vox.*origin;
		M   = [vox(1) 0 0 off(1) ; 0 vox(2) 0 off(2) ; 0 0 vox(3) off(3) ; 0 0 0 1];
	else
		% Assume it is a MINC file
		V=spm_vol_minc(imagename);
		if ~isempty(V),
			M=V.mat;
		else,
			% Try Ecat format
			V=spm_vol_ecat7(imagename);
			if ~isempty(V),
				M=V.mat;
			else,
				error(['Can''t read matrix information from "' imagename '".']);
			end;
		end;
	end;
elseif (nargin == 2)
	% Delete the matrix if it exists
	spm_unlink(matname);
	M    = mat;
	if exist([spm_str_manip(imagename,'sd') '.hdr']) == 2,
		% vx   = sqrt(sum(M(1:3,1:3).^2));
		% orgn = round(-M(1:3,4)' ./ vx);
		% off  = -vx.*orgn;
		% mt   = [vx(1) 0 0 off(1) ; 0 vx(2) 0 off(2) ; 0 0 vx(3) off(3) ; 0 0 0 1];

		[dim vox scale typ offset origin descrip] = spm_hread([spm_str_manip(imagename,'sd') '.img']);
		if isempty(dim),
			error(sprintf('Can''t read header for "%s"\n', [spm_str_manip(imagename,'sd') '.img']));
		end
		if all(origin == 0), origin = (dim(1:3)+1)/2; end;
		if all(vox == 0), vox = [1 1 1]; end;
		off = -vox.*origin;
		mt  = [vox(1) 0 0 off(1) ; 0 vox(2) 0 off(2) ; 0 0 vox(3) off(3) ; 0 0 0 1];

		% only write the .mat file if necessary
		if (sum((mat(:) - mt(:)).*(mat(:) - mt(:))) > eps*eps*12)
			save(matname,'M','-v4');
		end
	else,
		% Write the file anyway
		save(matname,'M','-v4');
	end;
else
	error('Incorrect Usage.');
end

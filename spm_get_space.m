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
[pth,nam,ext] = fileparts(imagename);
t = find(ext==',');
n = 1;
if ~isempty(t),
	if length(t)==1,
		n1 = ext((t+1):end);
		if ~isempty(n1),
			n   = str2num(n1);
			ext = ext(1:(t-1));
		end;
	end;
end;
mfname = fullfile(pth,[nam '.mat']);

if nargin == 1,

	% If we can read info from a ".mat" file with the corresponding
	% name, then use this information.
	if exist(mfname) == 2,
		clear M
		str = load(mfname);
		if isfield(str,'mat') & n<=size(str.mat,3) & any(any(str.mat(:,:,n))),
			M = str.mat(:,:,n);
			return;
		elseif isfield(str,'M'),
			M = str.M;
			if spm_flip_analyze_images & (strcmp(ext,'.img') | strcmp(ext,'.hdr')),
				M = diag([-1 1 1 1])*M;
			end;
			return;
		end;
	end;

	M = mat_from_header(imagename);

elseif nargin == 2,
	%v = spm_vol(sprintf('%s,%d', imagename, n));
	v = spm_vol(imagename);
	if sum((mat(:) - v.mat(:)).*(mat(:) - v.mat(:))) > eps*eps*12,
		% only do something if matrices are not identical
		if exist(mfname)==2,
			clear Mo
			str = load(mfname);
			if isfield(str,'mat'),
				Mo = str.mat;
			elseif isfield(str,'M'),
				Mo = str.M;
				if spm_flip_analyze_images & (strcmp(ext,'.img') | strcmp(ext,'.hdr')),
					Mo = diag([-1 1 1 1])*Mo;
				end;
			end;
			Mo(:,:,n) = mat;
			mat       = fill_empty(Mo,imagename);
			M = mat(:,:,1);
			if spm_flip_analyze_images & (strcmp(ext,'.img') | strcmp(ext,'.hdr')),
				M = diag([-1 1 1 1])*M;
			end;
			try,
				save(mfname,'mat','M','-append');
			catch,
				save(mfname,'mat','M');
			end;
		else,
			clear Mo
			Mo(:,:,n) = mat;
			mat       = fill_empty(Mo,imagename);
			M = mat(:,:,1);
			if spm_flip_analyze_images  & (strcmp(ext,'.img') | strcmp(ext,'.hdr')),
				M = diag([-1 1 1 1])*M;
			end;
			save(mfname,'mat','M');
		end;
	end;
else,
	error('Incorrect Usage.');
end;
return;

function Mo = fill_empty(Mo,imagename)
todo = [];
for i=1:size(Mo,3),
	if ~any(any(Mo(:,:,i))),
		todo = [todo i];
	end;
end;
if ~isempty(todo),
	Mfill = mat_from_header(imagename);
	for i=1:length(todo),
		Mo(:,:,todo(i)) = Mfill;
	end;
end;
return;

function M = mat_from_header(imagename)
[pth,nam,ext] = fileparts(imagename);
hfname        = fullfile(pth,[nam '.hdr']);
if exist(hfname) == 2,
	% Read as much information as is stored in the ANALYZE header.
	% Assume transverse slices.
	hdr = spm_read_hdr(hfname);
	if isempty(hdr),
		error(sprintf('Can''t read header for "%s"\n', imagename));
	end
	if any(hdr.hist.origin(1:3)),
		origin = hdr.hist.origin(1:3);
	else,
		origin = (hdr.dime.dim(2:4)+1)/2;
	end;
	vox    = hdr.dime.pixdim(2:4);
	if all(vox == 0), vox = [1 1 1]; end;
	off    = -vox.*origin;
	M      = [vox(1) 0 0 off(1) ; 0 vox(2) 0 off(2) ; 0 0 vox(3) off(3) ; 0 0 0 1];
	if spm_flip_analyze_images, M = diag([-1 1 1 1])*M; end;

else
	% Assume it is a MINC file
	V = spm_vol_minc(imagename);
	if ~isempty(V),
		M = V.mat;
	else,
		% Try Ecat format
		V = spm_vol_ecat7(imagename);
		if ~isempty(V),
			M = V.mat;
		else,
			error(['Can''t read matrix information from "' imagename '".']);
		end;
	end;
end;
return;


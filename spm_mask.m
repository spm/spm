function spm_mask(P)
% Masks Images.
% FORMAT spm_mask(P)
% P     - matrix of input image filenames.
% The masked images are prepended with the prefix `m'.
%
% If any voxel in the series of images is zero or does
% not have a finite value, then that voxel is set to NaN
% or zero in all the images.
%
% Images sampled in different orientations and positions can
% be passed to the routine.  Providing the `.mat' files are
% correct, then these should be handled appropriately.
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin==0,
	P=spm_get(Inf,'.img','Images to mask');
end;

V=spm_vol(P);

M=V(1).mat;
m=prod(size(V));

% Create headers
VO=V;
for i=1:prod(size(VO)),
	p  = spm_str_manip(VO(i).fname, 'd');
	q  = max([find(p == '/') 0]);
	q  = [p(1:q) 'm' p((q + 1):length(p))];
	VO(i).fname    = q;

	VO(i).descrip  = 'Masked';
	VO(i).mat      = VO(1).mat;
	VO(i).dim(1:3) = VO(1).dim(1:3);
	spm_create_image(VO(i));
end;

A = zeros([V(1).dim(1:2) m]);
spm_progress_bar('Init',V(1).dim(3),'Masking','planes completed')
for j=1:V(1).dim(3),

	% Load slice j from all images
	Mi=spm_matrix([0 0 j]);
	for i=1:m
		M1 = (M\V(i).mat)\Mi;
		A(:,:,i) = spm_slice_vol(V(i),M1,V(1).dim(1:2),[1 NaN]);
	end;

	% Get and apply the mask
	msk = get_mask(A);
	A   = apply_mask(A,msk);

	% Write the images.
	for i=1:m,
		tmp = A(:,:,i);
		spm_write_plane(VO(i),A(:,:,i),j);
	end;
	spm_progress_bar('Set',j);
end;
spm_progress_bar('Clear');
return;


function msk = get_mask(A)
	msk = find(any(~finite(A) | ~A,3));
return;

function A = apply_mask(A,msk)
	d = size(A);
	A = reshape(A,d(1)*d(2),d(3));
	for i=1:size(A,3),
		A(msk,:) = NaN;
	end;
	A = reshape(A,d);
return;

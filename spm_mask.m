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

% Extract voxel sizes and origin field from the matrix.
vx   = sqrt(sum(M(1:3,1:3).^2));
orig = M\[0 0 0 1]';
orig = round(orig(1:3)');

% Create headers and open files
for i=1:m,
	p  = spm_str_manip(V(i).fname, 'd');
	q  = max([find(p == '/') 0]);
	q  = [p(1:q) 'm' p((q + 1):length(p))];
	V(i).fp = fopen(q, 'w');
	V(i).q  = q;
	if V(i).fp == -1,
		open_error_message(q);
		error(['Error opening ' q '. Check that you have write permission.']);
	end;
	spm_hwrite(q,V(1).dim(1:3),vx,V(i).pinfo(1,1),V(i).dim(4),0,orig,'Masked');
end;

A=zeros([V(1).dim(1:2) m]);
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
		if any(V(i).dim(4) == [2 4 8]),
			msk = find(~finite(tmp));
			tmp(msk) = 0;
		end;
		l=fwrite(V(i).fp,tmp/V(i).pinfo(1),spm_type(V(i).dim(4)));
		if l~=prod(size(tmp)),
			spm_progress_bar('Clear');
			write_error_message(V(i).q);
			error(['Error writing ' V(i).q '. Check your disk space.']);
		end; 
	end;
	spm_progress_bar('Set',j);
end;
spm_progress_bar('Clear');

for i=1:m,
	fclose(V(i).fp);
end;
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


function open_error_message(q)
f=spm_figure('findwin','Graphics'); 
if ~isempty(f), 
	figure(f); 
	spm_figure('Clear','Graphics'); 
	spm_figure('Clear','Interactive'); 
	ax=axes('Visible','off','Parent',f); 
	text(0,0.60,'Error opening:', 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.55,spm_str_manip(q,'k40d'), 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.40,'  Please check that you have write permission.', 'FontSize', 16, 'Interpreter', 'none'); 
end
return

function write_error_message(q)
f=spm_figure('findwin','Graphics'); 
if ~isempty(f), 
	figure(f); 
	spm_figure('Clear','Graphics'); 
	spm_figure('Clear','Interactive'); 
	ax=axes('Visible','off','Parent',f); 
	text(0,0.60,'Error opening:', 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.55,spm_str_manip(q,'k40d'), 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.40,'  Please check that you have write permission.', 'FontSize', 16, 'Interpreter', 'none'); 
end
return;

function spm_mask(P1,P2, thresh)
% Masks Images.
% FORMAT spm_mask(P1,P2, thresh)
% P1     - matrix of input image filenames from which
%          to compute the mask.
% P2     - matrix of input image filenames on which
%          to apply the mask.
% thresh - optional threshold(s) for defining the mask.
% The masked images are prepended with the prefix `m'.
%
% If any voxel in the series of images is zero (for data types without
% a floating point representation) or does not have a finite value (for
% floating point and double precision images), then that voxel is set to
% NaN or zero in all the images.  If a threshold, or vector of
% thresholds is passed, then the masking is mased on voxels whos
% values are above all the thresholds.
%
% Images sampled in different orientations and positions can be passed
% to the routine.  Providing the `.mat' files are correct, then these
% should be handled appropriately.
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin==0,
	P1=spm_get(Inf,'.img','Images to compute mask from');
	P2=spm_get(Inf,'.img','Images to apply mask to');
end;
if nargin==1,
	P2 = P1;
end;
if nargin==4,
	if prod(size(thresh))==1,
		thresh = ones(prod(size(P1)),1)*thresh;
	else,
		if prod(size(P1)) ~= prod(size(thresh)),
			error('thresh vector is wrong size.');
		end;
		thresh=thresh(:);
	end;
end;

V1=spm_vol(P1);
V2=spm_vol(P2);

m1=prod(size(V1));
m2=prod(size(V2));

% Create headers
VO=V2;
for i=1:m2,
	[pth,nm,xt,vr] = fileparts(deblank(VO(i).fname));
	VO(i).fname    = fullfile(pth,['m' nm xt vr]);
	VO(i).descrip  = 'Masked';
	VO(i).mat      = VO(1).mat;
	VO(i).dim(1:3) = VO(1).dim(1:3);
	spm_create_image(VO(i));
end;
M   = VO(1).mat;
dim = VO(1).dim(1:3);

spm_progress_bar('Init',VO(1).dim(3),'Masking','planes completed')
for j=1:dim(3),

	msk = zeros(dim(1:2));
	Mi  = spm_matrix([0 0 j]);

	% Load slice j from all images
	for i=1:m1
		M1  = M\V1(i).mat\Mi;
		%if sum((M1(:)-Mi(:)).^2<eps) M1 = Mi; end;

		img = spm_slice_vol(V1(i),M1,dim(1:2),[0 NaN]);
		if nargin<3
			if ~spm_type(V1(i).dim(4),'nanrep'),
				msk = msk + (img~=0 & finite(img));
			else,
				msk = msk + finite(img);
			end;
		else,
			msk = msk + (img>=thresh & finite(img));
		end;
	end;

	msk = find(msk~=m1);

	% Write the images.
	for i=1:m2,
		M1  = M\V2(i).mat\Mi;
		img = spm_slice_vol(V2(i),M1,dim(1:2),[1 0]);
		img(msk) = NaN;
		spm_write_plane(VO(i),img,j);
	end;

	spm_progress_bar('Set',j);
end;
spm_progress_bar('Clear');
return;

function spm_warp_ui
% Warp a pair of same subject images together
%
% Very little support can be provided for the warping routine, as it
% involves optimising a very nonlinear objective function.
% Also, don't ask what the best value for the regularisation is.
%_______________________________________________________________________
% %W% John Ashburner %E%

VG = spm_vol(spm_get(1, '*.img', 'Reference image'));
VF = spm_vol(spm_get(1, '*.img', 'Image to warp'));

reg = spm_input('Regularisation constant','+0', 'e', '4', 1);
nit = spm_input('Number of iterations','+0', 'n', '8', 1);

% Load the image pair as 8 bit.
%-----------------------------------------------------------------------
VG.uint8 = loaduint8(VG); % Template
VF.uint8 = loaduint8(VF); % Surce image

% Try loading pre-existing deformation fields.  Otherwise, create
% deformation fields from uniform affine transformations.
%-----------------------------------------------------------------------
P = str2mat([prepend(VF.fname, 'y_') ',1'], ...
            [prepend(VF.fname, 'y_') ',2'], ...
            [prepend(VF.fname, 'y_') ',3']);
ok = 1;
if exist(prepend(VF.fname, 'y_'))~=2,
	ok = 0;
end;

if ok,
	VT = spm_vol(P);
	if any(abs(VT(1).mat\VG.mat - eye(4))>0.00001),
		error('Incompatible affine transformation matrices.');
	end;
	y1 = single(0);y1(VG.dim(1), VG.dim(2), VG.dim(3)) = 0;
	y2 = single(0);y2(VG.dim(1), VG.dim(2), VG.dim(3)) = 0;
	y3 = single(0);y3(VG.dim(1), VG.dim(2), VG.dim(3)) = 0;
	fprintf('Loading pre-existing deformation field\n');
	for i=1:VG.dim(3),
		% Note that this should read more like:
		% M = VT(1).mat\VG.mat*spm_matrix([0 0 i]);
		% Unfortunately, the warping algorithm does not cope with
		% NaNs.  Also, negative determinants can arise when a
		% deformation field is resampled using tri-linear
		% interpolation.

		M         = spm_matrix([0 0 i]);
		y1(:,:,i) = single(spm_slice_vol(VT(1),M,VG.dim(1:2),[1 NaN]));
		y2(:,:,i) = single(spm_slice_vol(VT(2),M,VG.dim(1:2),[1 NaN]));
		y3(:,:,i) = single(spm_slice_vol(VT(3),M,VG.dim(1:2),[1 NaN]));
	end;
	spm_affdef(y1,y2,y3,inv(VF.mat));
else,
	fprintf('Generating uniform affine transformation field\n');
	y1 = single(1:VG.dim(1))';
	y1 = y1(:,ones(VG.dim(2),1),ones(VG.dim(3),1));
	y2 = single(1:VG.dim(2));
	y2 = y2(ones(VG.dim(1),1),:,ones(VG.dim(3),1));
	y3 = reshape(single(1:VG.dim(3)),1,1,VG.dim(3));
	y3 = y3(ones(VG.dim(1),1),ones(VG.dim(2),1),:);
	spm_affdef(y1,y2,y3,VF.mat\VG.mat);
end;

% Voxel sizes
%-----------------------------------------------------------------------
vxg = sqrt(sum(VG.mat(1:3,1:3).^2))';if det(VG.mat(1:3,1:3))<0, vxg(1) = -vxg(1); end;
vxf = sqrt(sum(VF.mat(1:3,1:3).^2))';if det(VF.mat(1:3,1:3))<0, vxf(1) = -vxf(1); end;

% Do warping
%-----------------------------------------------------------------------
fprintf('Warping (iterations=%d regularisation=%g)\n', nit, reg);
spm_warp(VG.uint8,VF.uint8,y1,y2,y3,[vxg vxf],[nit,reg,1,1]);

% Convert mapping from voxels to mm
%-----------------------------------------------------------------------
spm_affdef(y1,y2,y3,VF.mat);

% Write the deformations
%-----------------------------------------------------------------------
VO         = VG;
VO.pinfo   = [1 0 0]';
VO.dim(4)  = spm_type('float');
VO.descrip = 'Deformation field';

VO.fname   = prepend(VF.fname, 'y_');
VO.n       = 1;
spm_write_vol(VO,y1);

VO.fname   = prepend(VF.fname, 'y_');
VO.n       = 2;
spm_write_vol(VO,y2);

VO.fname   = prepend(VF.fname, 'y_');
VO.n       = 3;
spm_write_vol(VO,y3);

return;
%_______________________________________________________________________

%_______________________________________________________________________
function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
if size(V.pinfo,2)==1 & V.pinfo(1) == 2,
	mx = 255*V.pinfo(1) + V.pinfo(2);
	mn = V.pinfo(2);
else,
	spm_progress_bar('Init',V.dim(3),...
		['Computing max/min of ' spm_str_manip(V.fname,'t')],...
		'Planes complete');
	mx = -Inf; mn =  Inf;
	for p=1:V.dim(3),
		img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
		mx  = max([max(img(:)) mx]);
		mn  = min([min(img(:)) mn]);
		spm_progress_bar('Set',p);
	end;
end;
spm_progress_bar('Init',V.dim(3),...
        ['Loading ' spm_str_manip(V.fname,'t')],...
        'Planes loaded');

udat = uint8(0);
udat(V.dim(1),V.dim(2),V.dim(3))=0;
for p=1:V.dim(3),
	img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
	udat(:,:,p) = uint8(round((img-mn)*((256-1)/(mx-mn))+1));
	spm_progress_bar('Set',p);
end;
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function out = prepend(in, pre)
[pth,nme,ext,ver] = fileparts(in);
out = fullfile(pth,[pre nme ext ver]);
return;
%_______________________________________________________________________

% Perform 3D stereotactic Normalization
% FORMAT spm_sn3d(P,matname,bb,Vox,params,spms)
% P         - image(s) to normalize
% matname   - name of file to store deformation definitions
% bb        - bounding box for normalized image
% Vox       - voxel sizes for normalized image
% params(1) - number of basis functions in X
% params(2) - "      "  "     "         "  Y
% params(3) - "      "  "     "         "  Z
% params(4) - number of iterations for elastic deformation
%             Setting any of these parameters to 0 will force
%             the program to perform only the affine
%             normalization.
% params(5) - smoothing for image (mm).
% params(6) - smoothness of deformation field.
% The routines attempt to match the input image to an
% optimum linear combination of template images. This
% provides additional flexibility in the type of input images
% which can be used. Typical template images consist of:
%     gray matter
%     white matter
%     scalp
%     striatum
%     ventricles
%     etc...
%
% First of all, a 12 parameter affine normalization stage
% is used to estimate the overall size, orientation etc of
% the image.
% Following this, an elastic deformation is computed which
% will match the image to the template.
% The deformation is in 3-dimensions, and is constrained
% to consist of a linear combination of basis functions.
% The basis functions chosen were the lowest frequency
% components of a 3 dimensional discrete cosine transform.
%
% The parameters for the affine transformation, and 3D
% basis function transformation are saved. The image
% is then resampled according these parameters, and
% written out (prefixing the filename with an "n").
%
% The code has been written, so that future versions should
% be able to handle multi-spectral input images.
%
% for a complete description of this approach see Friston et al (1994)
% ref: Friston et al (1994) The spatial registration and normalization of images
% HBM 0:00-00

% %W% John Ashburner MRCCU/FIL %E%

function spm_sn3d(P,matname,bb,Vox,params,spms)

% Map the template(s).
VG = [];
for i=1:size(spms,1)
	fname = spms(i,:);
	fname = fname(fname ~= ' ');
	VG = [VG spm_map(fname)];
	if (i==1)
		[dim,vox,scale,type,offset,origin,descrip] = spm_hread(fname);
	end
end

% Check for consistancy
if (size(VG,2)> 1 & any(diff(VG(1:6,:)')))
	error('Templates must have identical dimensions');
end

if (1==1)
% Map the image to normalize.
fprintf('smoothing..');
spm_smooth(P,'spm_sn3d_tmp.img',params(5));
VF = spm_map('spm_sn3d_tmp.img');
fprintf(' ..done\n');
else
VF = spm_map(P);
end

MF = spm_get_space(P(1,:));
MG = spm_get_space(spms(1,:));

if (1==1)
fprintf('Affine Normalization\n');
[Affine,scales] = spm_affine(VG,VF,MF\MG);
else
load(matname);
end

if (~any(params==0))
	fprintf('3D Cosine Transform Normalization\n');
	[Transform,Dims,scales] = spm_snbasis_map(VG,VF,Affine,params);
else
	Transform = [];
	Dims = [VG(1:3,1)' ; 0 0 0];
end

% Save parameters for future use.
Dims = [Dims ; vox ; origin ; VF(1:3,1)' ; VF(4:6,1)'];
mgc = 960209;
eval(['save ' matname ' mgc Affine Dims Transform MF MG']);

spm_unmap_vol(VF);

spm_write_sn('spm_sn3d_tmp.img',matname,bb,Vox,1);
delete spm_sn3d_tmp.img spm_sn3d_tmp.hdr spm_sn3d_tmp.mat
[tmp1,tmp2,tmp3,tmp4,tmp5,origin2,tmp6] = spm_hread('nspm_sn3d_tmp.img');

% Map the normalized volumes
VN = spm_map('nspm_sn3d_tmp.img');

% Do the display stuff
spm_clf;
centre = [10 0 0];

VGT = VG;
VGT(7,:) = VGT(7,:).*scales(:)';
figure(3); spm_clf
spm_orthviews(VGT,centre,origin ,bb,1,[0. 0. 1. .5],'Template');
spm_orthviews(VN ,centre,origin2,bb,1,[0. .5 1. .5],'Normalized');
drawnow;
delete nspm_sn3d_tmp.img nspm_sn3d_tmp.hdr nspm_sn3d_tmp.mat
spm_print

for v=[VG VN]
	spm_unmap_vol(v);
end

fprintf('Done\n');

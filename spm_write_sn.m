function spm_write_sn(P,matname,bb,Vox,Hold)
% Write Out Normalized Images.
% FORMAT spm_write_sn(P,matname,bb,Vox, Hold)
% P         - Image to transform.
% matname   - File containing transformation information.
% bb        - Bounding box (mm).
% Vox       - Required voxel size (mm).
% Hold      - Sampling method (see spm_sample_vol).
%_______________________________________________________________________
% %W% John Ashburner MRCCU/FIL %E%

[Dims,Affine,MF,MG,Tr] = load_params(matname);
V = spm_vol(P);

% The old voxel size and origin notation is used here.
% This requires that the position and orientation
% of the template is transverse.  It would not be
% straitforward to account for templates that are
% in different orientations because the basis functions
% would no longer be seperable.  The seperable basis
% functions mean that computing the deformation field
% from the parameters is much faster.

x = (bb(1,1):Vox(1):bb(2,1))/Dims(3,1) + Dims(4,1);
y = (bb(1,2):Vox(2):bb(2,2))/Dims(3,2) + Dims(4,2);
z = (bb(1,3):Vox(3):bb(2,3))/Dims(3,3) + Dims(4,3);

Dim = [length(x) length(y) length(z)];

X = x'*ones(1,Dim(2));
Y = ones(Dim(1),1)*y;

if (prod(size(Tr)) == 0),
	affine_only = 1;
	basX = 0; tx = 0;
	basY = 0; ty = 0;
	basZ = 0; tz = 0;
else,
	affine_only = 0;
	basX = spm_dctmtx(Dims(1,1),size(Tr,1),x-1);
	basY = spm_dctmtx(Dims(1,2),size(Tr,2),y-1);
	basZ = spm_dctmtx(Dims(1,3),size(Tr,3),z-1);
end;


spm_progress_bar('Init',length(z),'Computing available voxels','planes completed');
msk = cell(length(z),1);
for j=1:length(z),   % Cycle over planes
	Count = zeros(Dim(1),Dim(2));
	if affine_only,
		% Generate a mask for where there is data for all images
		%----------------------------------------------------------------------------
		for i=1:prod(size(V)),
			[X2,Y2,Z2] = mmult(X,Y,z(j),V(i).mat\MF*Affine);
			Count      = Count + getmask(X2,Y2,Z2,V(i).dim(1:3));
		end;
	else,
		% Nonlinear deformations
		%----------------------------------------------------------------------------
		tx = reshape(reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
		ty = reshape(reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
		tz = reshape(reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
		X1 = X    + basX*tx*basY';
		Y1 = Y    + basX*ty*basY';
		Z1 = z(j) + basX*tz*basY';

		% Generate a mask for where there is data for all images
		%----------------------------------------------------------------------------
		for i=1:prod(size(V)),
			[X2,Y2,Z2] = mmult(X1,Y1,Z1,V(i).mat\MF*Affine);
			Count      = Count + getmask(X2,Y2,Z2,V(i).dim(1:3));
		end;
	end;
	msk{j} = find(Count ~= prod(size(V)));
	spm_progress_bar('Set',j);
end;

% Loop over volumes.  This is an attempt to improve the efficiency of resampling the
% data, since the memory management should be able to keep the whole volume paged
% into RAM.  The only disadvantage is that the nonlinear deformations need to be
% computed over and over again - but this is pretty fast.  In theory, the deformations
% could be stored for a whole volume, but this would require quite a lot of memory.

spm_progress_bar('Init',prod(size(V)),'Resampling','volumes completed');
VO=V;
for i=1:prod(size(V)),
	VO(i).fname    = prepend(V(i).fname,'n');
	origin         = round(-bb(1,:)./Vox + 1);
	off            = -Vox.*origin;
	VO(i).mat      = [Vox(1) 0 0 off(1) ; 0 Vox(2) 0 off(2) ; 0 0 Vox(3) off(3) ; 0 0 0 1];
	VO(i).dim(1:3) = Dim;
	VO(i).descrip  = ['spm - 3D normalized'];
	spm_create_image(VO(i));

	if affine_only,
		for j=1:length(z),   % Cycle over planes
			[X2,Y2,Z2] = mmult(X,Y,z(j),V(i).mat\MF*Affine);
			tmp = spm_sample_vol(V(i),X2,Y2,Z2,Hold);
			tmp(msk{j}) = NaN;
			spm_write_plane(VO(i),tmp,j);
			if prod(size(V))<5, spm_progress_bar('Set',i-1+j/length(z)); end;
		end;
	else, % if ~affine_only,
		for j=1:length(z),   % Cycle over planes
			% Nonlinear deformations
			%----------------------------------------------------------------------------
			tx = reshape(reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			ty = reshape(reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			tz = reshape(reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3)) *basZ(j,:)', size(Tr,1), size(Tr,2) );
			X1 = X    + basX*tx*basY';
			Y1 = Y    + basX*ty*basY';
			Z1 = z(j) + basX*tz*basY';

			[X2,Y2,Z2] = mmult(X1,Y1,Z1,V(i).mat\MF*Affine);
			tmp = spm_sample_vol(V(i),X2,Y2,Z2,Hold);
			tmp(msk{j}) = NaN;
			spm_write_plane(VO(i),tmp,j);
			if prod(size(V))<5, spm_progress_bar('Set',i-1+j/length(z)); end;
		end;
	end;
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [Dims,Affine,MF,MG,Tr] = load_params(matname)
load(deblank(matname))
if (exist('mgc') ~= 1)
	error(['Matrix file ' matname ' is the wrong type.']);
end
if (mgc ~= 960209)
	error(['Matrix file ' matname ' is the wrong type.']);
end

% For evaluation using affine component only
%----------------------------------------------------------------------------
if 0,
	disp('Only using affine component');
	Dims(2,:)=[0 0 0];
	Tr=[];
end;
Tr = reshape(Transform,[Dims(2,:) 3]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function q = prepend(p,n)
p  = spm_str_manip(p, 'd');
q  = max([find(p == spm_platform('sepchar')) 0]);
q  = [p(1:q) n p((q + 1):length(p))];
return;
%_______________________________________________________________________

%_______________________________________________________________________
function Mask = getmask(X2,Y2,Z2,dim)
% Find range of slice
tiny = 5e-2;
Mask =        (X2 >= (1-tiny) & X2 <= (dim(1)+tiny));
Mask = Mask & (Y2 >= (1-tiny) & Y2 <= (dim(2)+tiny));
Mask = Mask & (Z2 >= (1-tiny) & Z2 <= (dim(3)+tiny));
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [X2,Y2,Z2] = mmult(X1,Y1,Z1,Mult);
if length(Z1) == 1,
	X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + (Mult(1,3)*Z1 + Mult(1,4));
	Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + (Mult(2,3)*Z1 + Mult(2,4));
	Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + (Mult(3,3)*Z1 + Mult(3,4));
else,
	X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
	Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
	Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
end;
return;
%_______________________________________________________________________

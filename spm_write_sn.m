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
	Transform=[];
end;


V = spm_vol(P);

% MP         maps this_image -> original_mm
% MG         maps template   -> normalised_mm
% MF         maps orig_image -> original_mm
% Transform  maps template   -> orig_image
%
% We want     template   -> this_image
%
%       Transform          MF          inv(MP)
% template -> orig_image -> original_mm -> this_image
%

% from mm space to voxel space.
x = (bb(1,1):Vox(1):bb(2,1))/Dims(3,1) + Dims(4,1);
y = (bb(1,2):Vox(2):bb(2,2))/Dims(3,2) + Dims(4,2);
z = (bb(1,3):Vox(3):bb(2,3))/Dims(3,3) + Dims(4,3);

Dim = [length(x) length(y) length(z)];

X = x'*ones(1,Dim(2));
Y = ones(Dim(1),1)*y;

if (prod(Dims(2,:)) == 0),
	affine_only = 1;
	basX = 0; tx = 0;
	basY = 0; ty = 0;
	basZ = 0; tz = 0;
else
	affine_only = 0;
	basX = spm_dctmtx(Dims(1,1),Dims(2,1),x-1);
	basY = spm_dctmtx(Dims(1,2),Dims(2,2),y-1);
	basZ = spm_dctmtx(Dims(1,3),Dims(2,3),z-1);
end

VO=V;
for i=1:prod(size(V)),
	% prefix filenames with `n'.
	p  = spm_str_manip(V(i).fname, 'd');
	q  = max([find(p == '/') 0]);
	q  = [p(1:q) 'n' p((q + 1):length(p))];

	VO(i).fname    = q;
	origin         = round(-bb(1,:)./Vox + 1);
	off            = -Vox.*origin;
	VO(i).mat      = [Vox(1) 0 0 off(1) ; 0 Vox(2) 0 off(2) ; 0 0 Vox(3) off(3) ; 0 0 0 1];
	VO(i).dim(1:3) = Dim;
	VO(i).descrip  = ['spm - 3D normalized'];
	spm_create_image(VO(i));
end;

% Start progress plot
%----------------------------------------------------------------------------
spm_progress_bar('Init',length(z),'Resampling','planes completed');

% Cycle over planes
%----------------------------------------------------------------------------
for j=1:length(z),
	% Nonlinear deformations
	%----------------------------------------------------------------------------
	if (~affine_only)
		% 2D transforms for each plane
		tx = reshape( reshape(Transform(:,1),Dims(2,1)*Dims(2,2),Dims(2,3)) *basZ(j,:)', Dims(2,1), Dims(2,2) );
		ty = reshape( reshape(Transform(:,2),Dims(2,1)*Dims(2,2),Dims(2,3)) *basZ(j,:)', Dims(2,1), Dims(2,2) );
		tz = reshape( reshape(Transform(:,3),Dims(2,1)*Dims(2,2),Dims(2,3)) *basZ(j,:)', Dims(2,1), Dims(2,2) );

		X1 = X    + basX*tx*basY';
		Y1 = Y    + basX*ty*basY';
		Z1 = z(j) + basX*tz*basY';
	end


	% Generate a mask for where there is data for all images
	%----------------------------------------------------------------------------
	Count = zeros(Dim(1),Dim(2));
	for i=1:size(P,1)
		Mult = V(i).mat\MF*Affine;
		if (~affine_only)
			X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
			Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
			Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
		else
			X2= Mult(1,1)*X + Mult(1,2)*Y + (Mult(1,3)*z(j) + Mult(1,4));
			Y2= Mult(2,1)*X + Mult(2,2)*Y + (Mult(2,3)*z(j) + Mult(2,4));
			Z2= Mult(3,1)*X + Mult(3,2)*Y + (Mult(3,3)*z(j) + Mult(3,4));
		end

		tiny = 5e-2;
		% Find range of slice

		Mask =        (X2 >= (1-tiny) & X2 <= (V(i).dim(1)+tiny));
		Mask = Mask & (Y2 >= (1-tiny) & Y2 <= (V(i).dim(2)+tiny));
		Mask = Mask & (Z2 >= (1-tiny) & Z2 <= (V(i).dim(3)+tiny));

		Count = Count + Mask;

	end
	Mask = (Count == size(P,1));

	for i=1:size(P,1)
		% Sample each volume
		%----------------------------------------------------------------------------
		Mult = V(i).mat\MF*Affine;
		if (~affine_only)
			X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
			Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
			Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
		else
			X2= Mult(1,1)*X + Mult(1,2)*Y + (Mult(1,3)*z(j) + Mult(1,4));
			Y2= Mult(2,1)*X + Mult(2,2)*Y + (Mult(2,3)*z(j) + Mult(2,4));
			Z2= Mult(3,1)*X + Mult(3,2)*Y + (Mult(3,3)*z(j) + Mult(3,4));
		end
		d = spm_sample_vol(V(i),X2,Y2,Z2,Hold).*Mask; % Apply mask
		spm_write_plane(VO(i),d,j);
	end
	spm_progress_bar('Set',j);
end
spm_progress_bar('Clear');
return;

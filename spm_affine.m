% Affine Normalization
% FORMAT [Affine,remainder] = spm_affine(VG,VF,Estimate);
% VG        - Vector of memory mapped templates.
% VF        - Memory mapped image to normalize.
% Estimate  - Starting estimate (optional).
% Affine    - The 4x4 transformation (in voxel space).
% remainder - The scaling for each of the templates required
%             to get the best match to the normalized image.
%
% spm_affine performs a 12 parameter normalization, using
% translations, rotations, zooms, and skew.
% The images should all have roughly an 8mm FWHM smoothness.

% %W% John Ashburner MRCCU/FIL %E%


function [Affine,remainder] = spm_affine(VG,VF,Estimate)

if (nargin ~= 3)
	% Starting estimate
	z = VG(4:6,1)./VF(4:6,1);
	Estimate = diag([z; 1]);
	Estimate(:,4) = [(VF(1:3,1)-VG(1:3,1).*z)/2; 1];
end

Affine = Estimate;

% X - coordianates of templates
X = sparse(1:VG(1), 1:VG(1), 1:VG(1), VG(1), VG(1)) * ones(VG(1), VG(2));

% Y - coordinates of templates
Y = ones(VG(1), VG(2)) * sparse(1:VG(2), 1:VG(2), 1:VG(2), VG(2), VG(2));

% Sample about every 4 mm
samp=4;
skipx = max([round(samp/VG(4)) 1]);
skipy = max([round(samp/VG(5)) 1]);
skipz = max([round(samp/VG(6)) 1]);
mask0 = find((rem(X,skipx) == 0) & (rem(Y,skipy) == 0));

% Coordinates of corners of template.
corners = [
1     1     1     1
1     1     VG(1) 1
1     VG(2) 1     1
1     VG(2) VG(1) 1
VG(3) 1     1     1
VG(3) 1     VG(1) 1
VG(3) VG(2) 1     1
VG(3) VG(2) VG(1) 1]';

% Coordinates to which corners of template transform.
corners1 = Affine*corners;
corners1 = diag(VF(4:6,1))*corners1(1:3,:);

np = 12;
scales = ones(size(VG,2),1);

for iter=1:64
	fprintf('iteration # %d:', iter);

	alpha = zeros(np+size(VG,2),np+size(VG,2));
	beta  = zeros(np+size(VG,2),1);

	Mat = Affine(1:3,1:3)';

	for p=1:skipz:VG(3),

		XM = X(mask0);
		YM = Y(mask0);

		% Transformed template coordinates.
		X1= Affine(1,1)*XM + Affine(1,2)*YM + (Affine(1,3)*p + Affine(1,4));
		Y1= Affine(2,1)*XM + Affine(2,2)*YM + (Affine(2,3)*p + Affine(2,4));
		Z1= Affine(3,1)*XM + Affine(3,2)*YM + (Affine(3,3)*p + Affine(3,4));

		% Only resample from within the volume VF.
		mask1 = find((Z1>=1) & (Z1<VF(3)-0.01) ...
			   & (Y1>=1) & (Y1<VF(2)-0.01) ...
			   & (X1>=1) & (X1<VF(1)-0.01));

		% Don't waste time on an empty plane.
		if (length(mask1>0))

			% Only resample from within the volume VF.
			if (length(mask1) ~= length(mask0))
				X1 = X1(mask1);
				Y1 = Y1(mask1);
				Z1 = Z1(mask1);
				XM = XM(mask1);
				YM = YM(mask1);
			end
			ZM = zeros(size(mask1))+p;

			dFdM = zeros(size(mask1,1),np+size(VG,2));

			% Sample image to normalise & get local NEGATIVE derivatives
			F  =     spm_sample_vol(VF, X1    , Y1    , Z1    , 1 );

			% Sample referance image(s)
			for i=1:size(VG,2)
				dFdM(:,np+i)=spm_sample_vol(VG(:,i), XM, YM, ZM, 1)*scales(i) - F;
			end

			% Not yet 100% sure about this.
			dx =(F - spm_sample_vol(VF, X1+.01, Y1    , Z1    , 1 ))/(.01);
			dy =(F - spm_sample_vol(VF, X1    , Y1+.01, Z1    , 1 ))/(.01);
			dz =(F - spm_sample_vol(VF, X1    , Y1    , Z1+.01, 1 ))/(.01);

			% Generate Design Matrix
			dFdM(:,1:np) = [ XM.*dx YM.*dx p*dx dx ...
					 XM.*dy YM.*dy p*dy dy ...
					 XM.*dz YM.*dz p*dz dz ];


			% Most of the work
			alpha = alpha + spm_atranspa(dFdM);
			beta = beta + dFdM'*F;
			clear('dFdM');
		end
		fprintf('.');
	end
	% Least squares solution
	q = alpha\beta;
	Affine(1:3,1:4) = Affine(1:3,1:4) + reshape(q(1:np),4,3)';
	scales = scales + q((1:size(VG,2))+np);

	% If the maximum movement is less than 0.25 mm since the
	% last iteration - then finished.
	corners2 = corners1;
	corners1 = Affine*corners;
	corners1 = diag(VF(4:6,1))*corners1(1:3,:);
	movement = max(sqrt(sum(corners1-corners2).^2));
	fprintf(' %.3f\n', movement);

	if (movement < 0.25 & iter >= 12)
		break;
	end
end
remainder = q((1:size(VG,2)) + np);
fprintf('\n');


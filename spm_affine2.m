function params = spm_affine2(VG, VF, MG, MF, params, free, Hold, tol, samp)
% Affine Normalization/Coregistration.
% FORMAT params = spm_affine2(VG, VF, MG, MF, params, free, Hold, tol, samp)
% VG        - Vector of memory mapped targets.
% VF        - Memory mapped object images.
% MG        - Space of the targets
% MF        - Space of the objects.
%
%	Optional Parameters:
% params    - Starting estimates.
%             	defaults: [1 1 1 1 1 1 1 1 1 1 1 1] - fit all parameters.
% free      - Parameters to be fitted.
%             	defaults: [0 0 0 0 0 0 1 1 1 0 0 0].
% Hold      - Interpolation method.
%             	default: 1 - bilinear interpolation
% tol       - Tolerance (maximum shift is less than tol mm then terminate).
%             	default: 0.25 mm.
% samp      - Sample distance (sample template(s) roughly every samp mm).
%             	default: 4 mm.
%
%------------------------------------------------------------------------
% Simultaneously matches F1 to G1, F2 to G2 etc. The images are assumed
% to match without scaling any of them.
%
% Each iteration solves (for p1, p2 .. pn):
%
% / dF1/dp1 dF1/dp1 .. dF1/dpn \                    / (F1-G1) \
% | dF2/dp1 dF2/dp1 .. dF2/dpn | * [p1 p2 .. pn]' = | (F2-G2) |
% \ dFm/dp1 dFm/dp1 .. dFm/dpn /                    \ (Fm-Gm) /
%

% %W% John Ashburner FIL %E%

if (size(VF,2) ~= size(VG,2))
	error('There must be equal numbers of objects and targets');
end

% Set up default parameters.
%-----------------------------------------------------------------------
if (nargin<9)
	samp = 4;
	if (nargin<8)
		tol = 0.25;
		if (nargin<7)
			Hold = 1;
			if (nargin<6)
				free = [1 1 1 1 1 1 1 1 1 1 1 1];
				if (nargin<5)
					params = [0 0 0 0 0 0 1 1 1 0 0 0];
				end
			end
		end
	end
end

params = params(:);
free = find(free(:));

% coordianates of templates
%-----------------------------------------------------------------------
X = sparse(1:VG(1), 1:VG(1), 1:VG(1), VG(1), VG(1)) * ones(VG(1), VG(2));
Y = ones(VG(1), VG(2)) * sparse(1:VG(2), 1:VG(2), 1:VG(2), VG(2), VG(2));

% Sample about every samp mm
%-----------------------------------------------------------------------
skipx = max([round(samp/VG(4)) 1]);
skipy = max([round(samp/VG(5)) 1]);
skipz = max([round(samp/VG(6)) 1]);
mask0 = find((rem(X,skipx) == 0) & (rem(Y,skipy) == 0));

% Coordinates of corners of template.
%-----------------------------------------------------------------------
corners = [
1     1     1     1
1     1     VG(3) 1
1     VG(2) 1     1
1     VG(2) VG(3) 1
VG(1) 1     1     1
VG(1) 1     VG(3) 1
VG(1) VG(2) 1     1
VG(1) VG(2) VG(3) 1]';

% Coordinates to which corners of template transform.
%-----------------------------------------------------------------------
Mat = inv(inv(MG)*spm_matrix(params')*MF);
corners1 = Mat*corners;
corners1 = diag(VF(4:6,1))*corners1(1:3,:);

np = size(free,1);

tmp   = spm_matrix(params');
origd = det(tmp(1:3,1:3));

for iter=1:64
	fprintf('iteration # %d:', iter);

	alpha = zeros(12);
	beta  = zeros(12,1);

	Mat = inv(inv(MG)*spm_matrix(params')*MF);

	% rate of change of matrix elements with respect to parameters
	%-------------------------------------------------------------
	dMdP = zeros(12, np);
	t0  = Mat(1:3,1:4)';
	t0   = t0(:);
	for pp = 1:np;
		p = free(pp);
		tparams = params;
		tparams(p) = tparams(p)+0.00001;
		tmp = inv(inv(MG)*spm_matrix(tparams')*MF);
		tmp = tmp(1:3,1:4)';
		dMdP(:,p) = (tmp(:)-t0)/0.00001;
	end

	% Loop over planes
	%-------------------------------------------------------------
	for p=1:skipz:VG(3),

		XM = X(mask0);
		YM = Y(mask0);

		% Transformed template coordinates.
		%-----------------------------------------------------------------------
		X1= Mat(1,1)*XM + Mat(1,2)*YM + (Mat(1,3)*p + Mat(1,4));
		Y1= Mat(2,1)*XM + Mat(2,2)*YM + (Mat(2,3)*p + Mat(2,4));
		Z1= Mat(3,1)*XM + Mat(3,2)*YM + (Mat(3,3)*p + Mat(3,4));

		% Only resample from within the volume VF.
		%-----------------------------------------------------------------------
		mask1 = find((Z1>=1) & (Z1<VF(3)-0.01) ...
			   & (Y1>=1) & (Y1<VF(2)-0.01) ...
			   & (X1>=1) & (X1<VF(1)-0.01));


		if (length(mask1>0))	% Don't waste time on an empty plane.

			% Only resample from within the volume VF.
			%-----------------------------------------------------------------------
			if (length(mask1) ~= length(mask0))	
				X1 = X1(mask1);
				Y1 = Y1(mask1);
				Z1 = Z1(mask1);
				XM = XM(mask1);
				YM = YM(mask1);
			end
			ZM = zeros(size(mask1))+p;

			for j=1:size(VF,2)

				% Sample object image & get local NEGATIVE derivatives
				%-----------------------------------------------------------------------
				F  = spm_sample_vol(VF(:,j), X1 ,Y1, Z1, Hold);
				G  = spm_sample_vol(VG(:,j), XM, YM, ZM, Hold);

				dx =(F - spm_sample_vol(VF(:,j),X1+.001,Y1     ,Z1     ,Hold))/(.001);
				dy =(F - spm_sample_vol(VF(:,j),X1     ,Y1+.001,Z1     ,Hold))/(.001);
				dz =(F - spm_sample_vol(VF(:,j),X1     ,Y1     ,Z1+.001,Hold))/(.001);


				% Generate Design Matrix
				%-----------------------------------------------------------------------
				dFdM = [ XM.*dx YM.*dx p*dx dx ...
					 XM.*dy YM.*dy p*dy dy ...
					 XM.*dz YM.*dz p*dz dz ];

				% Update alpha'*alpha and alpha'*beta
				%-----------------------------------------------------------------------
				alpha = alpha + spm_atranspa(dFdM);
				beta  = beta  + dFdM'*(F-G);

				clear('dFdM');
			end
		end
		fprintf('.');
	end

	% Least squares solution. (works because: dF/dP = dF/dM * dM/dP)
	%-----------------------------------------------------------------------
	q = inv(dMdP'*alpha*dMdP)*(dMdP'*beta);
	params(free) = params(free) + q(1:np);

	Mat = inv(inv(MG)*spm_matrix(params')*MF);


	% A safety measure for the unlikely case that the orientation
	% of the images flips without anyone noticing.
	% (ie. one of the zooms becomes -ve).
	%-----------------------------------------------------------------------
	tmp  = spm_matrix(params');
	newd = det(tmp(1:3,1:3));
	if (sign(newd/origd) < 0)
		error('Coordinate system has flipped');
	end

	% If the maximum movement is less than tol mm since the
	% last iteration - then finished.
	%-----------------------------------------------------------------------
	corners2 = corners1;
	corners1 = Mat*corners;
	corners1 = diag(VF(4:6,1))*corners1(1:3,:);
	movement = max(sqrt(sum(corners1-corners2).^2));
	fprintf(' %.3f\n', movement);
	if (movement < tol)
		fprintf('\n');
		return;
	end
end

error('This registration doesn''t look like it''s going to work');

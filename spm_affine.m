function [params,scales] = spm_affine(VG,VF,MG,MF,params,free,Hold,tol,samp)
% Affine Normalization
% FORMAT [params,scales] = spm_affine(VG,VF,MG,MF,params,free,Hold, tol, samp)
% VG        - Vector of memory mapped templates.
% VF        - Memory mapped image to normalize.
% MG        - Space of the templates
% MF        - Space of the image to normalise.
%
%	Optional Parameters:
% params    - Starting estimates.
%             	defaults: [0 0 0 0 0 0 1 1 1 0 0 0].
% free      - Parameters to be fitted.
%             	defaults: [1 1 1 1 1 1 1 1 1 1 1 1] - fit all parameters.
% Hold      - Interpolation method.
%             	default: 1 - bilinear interpolation
% tol       - Tolerance (maximum shift is less than tol mm then terminate).
%             	default: 0.4 mm.
% samp      - Sample distance (sample template(s) roughly every samp mm.
%             	default: 4 mm.
%
%-----------------------------------------------------------------------
%
% Determine the affine transformation which best fits image F to a linear
% combination of G1, G2, etc..
% Each iteration solves (for p1, p2 .. pn):
% [ G1 G2 .. dF/dp1 dF/dp1 ..] * [s1 s2 .. p1 p2 ..]' = [F]
%
%
% The transformation from one space to another can be obtained by:
% 	MG\spm_matrix(params')*MF

% %W% John Ashburner FIL %E%


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

% Coordianates of templates
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

np = prod(size(params));
np = size(free,1);

tmp   = spm_matrix(params');
origd = det(tmp(1:3,1:3));

fprintf('Determining Affine Transformation\n');
countdown = 0;
iter = 1;
while (iter <= 64 & countdown <= 4)

	fprintf('iteration # %d:', iter);

	alpha = zeros(12+size(VG,2),12+size(VG,2));
	beta  = zeros(12+size(VG,2),1);

	Mat = inv(inv(MG)*spm_matrix(params')*MF);

	% rate of change of matrix elements with respect to parameters
	%-----------------------------------------------------------------------
	dMdP = zeros(12+size(VG,2),(np+size(VG,2)));
	tmp = Mat(1:3,1:4)';
	t0 = [tmp(:); zeros(size(VG,2),1)];
	for pp = 1:np;
		p = free(pp);
		tparams = params;
		tparams(p) = tparams(p)+0.01;
		tmp = inv(inv(MG)*spm_matrix(tparams')*MF);
		tmp = tmp(1:3,1:4)';
		dMdP(:,p) = ([tmp(:); zeros(size(VG,2),1)]-t0)/0.01;
	end
	dMdP(:,(np+1):(np+size(VG,2))) = [zeros(12,size(VG,2)); eye(size(VG,2))];

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

		% Don't waste time on an empty plane.
		%-----------------------------------------------------------------------
		if (length(mask1>0))

			% Only resample from within the volume VF.
			%-----------------------------------------------------------------------
			if (length(mask1) ~= length(mask0))
				X1 = X1(mask1);
				Y1 = Y1(mask1);
				Z1 = Z1(mask1);
				XM = XM(mask1);
				YM = YM(mask1);
			end

			dFdM = zeros(size(mask1,1),12+size(VG,2));

			% Sample image to normalise & get local NEGATIVE derivatives
			%-----------------------------------------------------------------------
			F  =      spm_sample_vol(VF,X1    ,Y1    ,Z1    ,Hold);
			dx = (F - spm_sample_vol(VF,X1+.01,Y1    ,Z1    ,Hold))/(.01);
			dy = (F - spm_sample_vol(VF,X1    ,Y1+.01,Z1    ,Hold))/(.01);
			dz = (F - spm_sample_vol(VF,X1    ,Y1    ,Z1+.01,Hold))/(.01);

			% Generate Design Matrix
			%-----------------------------------------------------------------------
			dFdM(:,1:12) = [ XM.*dx YM.*dx p*dx dx ...
					 XM.*dy YM.*dy p*dy dy ...
					 XM.*dz YM.*dz p*dz dz ];

			% Sample referance image(s)
			%-----------------------------------------------------------------------
			ZM = zeros(size(mask1))+p;
			for i=1:size(VG,2)
				dFdM(:,12+i)=spm_sample_vol(VG(:,i), XM, YM, ZM, Hold);
			end


			% Most of the work
			%-----------------------------------------------------------------------
			alpha = alpha + spm_atranspa(dFdM);
			beta = beta + dFdM'*F;
			clear('dFdM');
		end
		fprintf('.');
	end

	% Least squares solution
	%-----------------------------------------------------------------------
	q = pinv(dMdP'*alpha*dMdP)*(dMdP'*beta);
	params(free) = params(free) + q(1:np);

	% A safety measure for the unlikely case that the orientation
	% of the images flips without anyone noticing.
	%-----------------------------------------------------------------------
	tmpM  = spm_matrix(params');
	newd = det(tmpM(1:3,1:3));
	if (sign(newd/origd) < 0)
		error('Coordinate system has flipped');
	end


	% If the maximum movement is less than tol mm since the
	% last iteration - then repeat for another few iterations.
	%-----------------------------------------------------------------------
	Mat = inv(inv(MG)*tmpM*MF);
	corners2 = corners1;
	corners1 = Mat*corners;
	corners1 = diag(VF(4:6,1))*corners1(1:3,:);
	movement = max(sqrt(sum(corners1-corners2).^2));
	fprintf(' %.3f\n', movement);

	if (movement < tol)
		countdown = countdown + 1;
	else
		countdown = 0;
	end

	iter = iter + 1;
end

scales = q((1:size(VG,2))+np);
fprintf('\n');


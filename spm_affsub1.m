function [alpha, beta, chi2] = spm_affsub1(VG,VF,MG,MF,Hold,samp,P)
% Generate A'*A and A'*b and \Chi^2 for affine image registration.
% FORMAT [alpha, beta, chi2] = spm_affsub1(VG,VF,MG,MF,P,Hold,samp)
% VG        - Vector of memory mapped template image(s).
% VF        - Memory mapped object image.
% MG        - Space of the template image(s).
% MF        - Space of the object image.
% Hold      - Interpolation method.
% samp      - frequency (in mm) of sampling.
% P         - Current parameter estimates.
% 
% alpha     - A'*A
% beta      - A'*b
% chi2      - Residual sum of squares divided by number of voxels sampled.
%__________________________________________________________________________
%
% Compare this subroutine with "mrqcof" from "Numerical Recipes".
% The parameters are:
%    P(1)  - x translation
%    P(2)  - y translation
%    P(3)  - z translation
%    P(4)  - x rotation about - {pitch} (radians)
%    P(5)  - y rotation about - {roll}  (radians)
%    P(6)  - z rotation about - {yaw}   (radians)
%    P(7)  - x scaling
%    P(8)  - y scaling
%    P(9)  - z scaling
%    P(10) - x affine
%    P(11) - y affine
%    P(12) - z affine
%    P(13) - scale required for image G(1) to best fit image F.
%
% With more than 13 parameters, then parameters 13 onwards describe a
% linear combination of the template images.
%
%__________________________________________________________________________
% %W% John Ashburner FIL %E%

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


alpha = zeros(12+size(VG,2),12+size(VG,2));
beta  = zeros(12+size(VG,2),1);

Mat = inv(inv(MG)*spm_matrix(P')*MF);

% rate of change of matrix elements with respect to parameters
%-----------------------------------------------------------------------
dMdP = zeros(12+size(VG,2),(12+size(VG,2)));
tmp = Mat(1:3,1:4)';
t0 = [tmp(:); zeros(size(VG,2),1)];
for pp = 1:12;
	tP = P;
	tP(pp) = tP(pp)+0.01;
	tmp = inv(inv(MG)*spm_matrix(tP')*MF);
	tmp = tmp(1:3,1:4)';
	dMdP(:,pp) = ([tmp(:); zeros(size(VG,2),1)]-t0)/0.01;
end
dMdP(:,(1:size(VG,2))+12) = [zeros(12,size(VG,2)); eye(size(VG,2))];



chi2 = 0;
n    = 0;

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
			tmp          = spm_sample_vol(VG(:,i), XM, YM, ZM, Hold);
			dFdM(:,12+i) = tmp;
			F            = F - tmp*P(i+12);
		end

		chi2 = chi2 + sum(sum(F.*F));
		n    = n    + prod(size(F));

		% Most of the work
		%-----------------------------------------------------------------------
		alpha = alpha + spm_atranspa(dFdM);
		beta = beta + dFdM'*F;
		clear('dFdM');
	end
	%fprintf('.');
end
%fprintf('\n');

alpha  = dMdP'*alpha*dMdP;
beta   = dMdP'*beta;
chi2   = chi2/(n - 12 - size(VG,2));

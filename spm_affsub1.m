function [alpha, beta, chi2, W] = spm_affsub1(VG,VF,MG,MF,Hold,samp,P,flg,minW)
% Generate A'*A and A'*b and \Chi^2 for affine image registration.
% FORMAT [alpha, beta, chi2, W] = spm_affsub1(VG,VF,MG,MF,Hold,samp,P,minW)
% VG        - Vector of memory mapped template image(s).
% VF        - Memory mapped object image.
% MG        - Space of the template image(s).
% MF        - Space of the object image.
% Hold      - Interpolation method.
% samp      - frequency (in mm) of sampling.
% P         - Current parameter estimates.
% minW      - previous minimum smoothness estimate.
% 
% alpha     - A'*A
% beta      - A'*b
% chi2      - Residual sum of squares.
% flg       - flag - a value of one means use the derivatives of F and G
%             rather than just F.
% W         - smoothness estimate.
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
% Parameters 13 onwards describe a linear combination of the
% template images.
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

% Convert parameters to affine transformation matrix
%-----------------------------------------------------------------------
Mat = inv(inv(MG)*spm_matrix(P')*MF);

% rate of change of matrix elements with respect to parameters
%-----------------------------------------------------------------------
dMdP = zeros(12+size(VG,2),(12+size(VG,2)));
tmp = Mat(1:3,1:4)';
t0 = [tmp(:); zeros(size(VG,2),1)];
for pp = 1:12;
	tP = P;
	tP(pp) = tP(pp)+0.001;
	tmp = inv(inv(MG)*spm_matrix(tP')*MF);
	tmp = tmp(1:3,1:4)';
	dMdP(:,pp) = ([tmp(:); zeros(size(VG,2),1)]-t0)/0.001;
end
dMdP(:,(1:size(VG,2))+12) = [zeros(12,size(VG,2)); eye(size(VG,2))];

% Initialise variables
%-----------------------------------------------------------------------
alpha = zeros(12+size(VG,2),12+size(VG,2));
beta  = zeros(12+size(VG,2),1);
chi2  = 0;
dch2  = [0 0 0];
n     = 0;

for p=1:skipz:VG(3),	% loop over planes
	XM = X(mask0);
	YM = Y(mask0);

	% Transformed template coordinates.
	%-----------------------------------------------------------------------
	X1= Mat(1,1)*XM + Mat(1,2)*YM + (Mat(1,3)*p + Mat(1,4));
	Y1= Mat(2,1)*XM + Mat(2,2)*YM + (Mat(2,3)*p + Mat(2,4));
	Z1= Mat(3,1)*XM + Mat(3,2)*YM + (Mat(3,3)*p + Mat(3,4));

	% Only resample from within the volume VF.
	%-----------------------------------------------------------------------
	t = 4.9e-2;
	mask1 = find((Z1>=1-t) & (Z1<=VF(3)+t) ...
		   & (Y1>=1-t) & (Y1<=VF(2)+t) ...
		   & (X1>=1-t) & (X1<=VF(1)+t));

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
		ZM = zeros(size(mask1))+p;

		% Rate of change of residuals w.r.t parameters
		%-----------------------------------------------------------------------
		dResdM = zeros(size(mask1,1),12+size(VG,2));

		% Sample object image & get local derivatives
		%-----------------------------------------------------------------------
		[F,dxF,dyF,dzF] = spm_sample_vol(VF, X1, Y1, Z1, Hold);

		% Sample referance image(s) and derivatives
		%-----------------------------------------------------------------------
		for i=1:size(VG,2)

			if flg==1 | nargout>=4,
				[Gi,dxt,dyt,dzt] = spm_sample_vol(VG(:,i), XM, YM, ZM, Hold);

				if i==1
					res = F   - Gi*P(i+12); % Residuals
					dxG = dxt*P(i+12);
					dyG = dyt*P(i+12);
					dzG = dzt*P(i+12);
				else
					res = res - Gi*P(i+12);
					dxG = dxG+dxt*P(i+12);
					dyG = dyG+dyt*P(i+12);
					dzG = dzG+dzt*P(i+12);
				end
			else,
				Gi = spm_sample_vol(VG(:,i), XM, YM, ZM, Hold);
				if i==1
					res = F   - Gi*P(i+12); % Residuals
				else
					res = res - Gi*P(i+12);
				end
			end;

			dResdM(:,12+i) = -Gi;

		end

		if (flg==1)
			% Use derivatives of F and G
			tmp = inv(Mat(1:3,1:3))';
			dG  = [dxG dyG dzG]*tmp';       % Rotate derivatives to space of F
			dx = 0.5*(dxF + dG(:,1));
			dy = 0.5*(dyF + dG(:,2));
			dz = 0.5*(dzF + dG(:,3));
		else
			% Do it the old way
			dx = dxF;
			dy = dyF;
			dz = dzF;
		end

		% Generate Design Matrix from rate of change of residuals wrt matrix
		% elements.
		%-----------------------------------------------------------------------
		dResdM(:,1:12) = [ XM.*dx YM.*dx p*dx dx ...
		                   XM.*dy YM.*dy p*dy dy ...
		                   XM.*dz YM.*dz p*dz dz ];

		% alpha = alpha + A'*A and beta = beta + A'*b
		%-----------------------------------------------------------------------
		alpha = alpha + spm_atranspa(dResdM);
		beta  = beta + dResdM'*res;
		clear dResdM

		% Assorted variables which are used later.
		%-----------------------------------------------------------------------
		chi2  = chi2 + res'*res;		% Sum of squares of residuals
		n     = n    + prod(size(F));		% Number of observations

		if nargout>=4,
			% Spatial derivatives of residuals derived from
			% (derivatives of F rotated to space of G) - (derivatives of G).
			%-----------------------------------------------------------------------
			tmp = Mat(1:3,1:3)';
			dF  = [dxF dyF dzF]*tmp';
			dx  = dF(:,1) - dxG;
			dy  = dF(:,2) - dyG;
			dz  = dF(:,3) - dzG;
			dch2  = dch2 + [dx'*dx dy'*dy dz'*dz];	% S.o.sq of derivs of residls
		end
	end
end

if n<=1
	f=spm_figure('findwin','Graphics');
	if ~isempty(f),
		figure(f);
		spm_figure('Clear','Graphics');
		spm_figure('Clear','Interactive');
		ax=axes('Visible','off','Parent',f);
		text(0,0.60,'There is not enough overlap in the', 'FontSize', 25, 'Interpreter', 'none');
		text(0,0.55,'    images to obtain a solution.', 'FontSize', 25, 'Interpreter', 'none');
		text(0,0.40,'  Please check that your header information is OK.','FontSize', 16, 'Interpreter', 'none');
	end

	error('There is not enough overlap of the images to obtain a solution');
end

% Smoothness estimate from residuals to determine number of
% independant observations.
%-----------------------------------------------------------------------
if nargout>=4,
	W      = (2*dch2/chi2).^(-.5).*VG(4:6,1)';
	msk    = find(~finite(W));
	W      = min([W;minW]);
	W(msk) = 1;
	skips  = [skipx skipy skipz].*VG(4:6,1)';		% sample distances (mm)
	skips(msk)  = 1;

%fprintf('fwhm = [%g %g %g]\n', W*sqrt(8*log(2)));
%df    = (n - size(beta,1))*prod(erf(2^(-3/2) * skips./W(msk)));	% possible alternative

	smo    = prod(min(skips./(W*sqrt(2*pi)),[1 1 1])); % fmri revisited
else
	smo = 1;
end

df     = (n - size(beta,1))*smo;

%fprintf('%g\n',chi2/df);

% Compute alpha and beta
%-----------------------------------------------------------------------
chi2  = chi2/df;
alpha = dMdP'*alpha*dMdP/chi2;
beta  = dMdP'*beta      /chi2;

return;

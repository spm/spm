function params = spm_affsub3(mode, VG, VF, Hold, samp, params,VW,VW2)
% Highest level subroutine involved in affine transformations.
% FORMAT params = spm_affsub3(mode, VG, VF, Hold, samp, params,VW,VW2)
%
% mode      - Mode of action.
% VG        - Handles of template images (see spm_vol).
% VF        - Handles of object images.
% Hold      - Interpolation method.
% samp      - Frequency (in mm) of sampling.
% params    - Parameter estimates.
%
% optional:
% VW        - Handle of weight image.
% VW2       - Handle of weight image for object image(s)
%__________________________________________________________________________
%
% Currently mode must be one of the following:
%	'register1'
%		This is for use in Multimodal coregistration.
%		Each F is mapped to one G (with scaling), but the
%		rigid body components differ between the two sets
%		of registrations.
%	'rigid1'
%		Rigid body registration.
%		Each F is mapped to one G, without scaling.
%	'rigid2'
%		Rigid body registration.
%		Each F is mapped to one G, with scaling.
%	'rigid3'
%		Rigid body registration.
%		Each F is mapped to a linear combination of Gs.
%	'affine1'
%		Affine normalisation.
%		Each F is mapped to one G, without scaling.
%	'affine2'
%		Affine normalisation.
%		Each F is mapped to one G, with scaling.
%	'affine3'
%		Affine normalisation.
%		Each F is mapped to a linear combination of Gs.
%
%	'2d1'
%		For 2d rigid-body registration.
%__________________________________________________________________________
% %W% John Ashburner FIL %E%

if nargin<5 | nargin>8,
	error('Incorrect usage.');
end;

global sptl_Ornt

% Covariance matrix of affine normalization parameters. Mostly assumed
% to be diagonal, but covariances between the three zooms is also
% accounted for.
% Data determined from 51 normal brains.
%-----------------------------------------------------------------------
c11 = eye(3)*10000;                % Allow stdev of 100mm in all directions.
c22 = eye(3)*0.046;                % 7 degrees standard deviations.
c33 = [0.0021 0.0009 0.0013        % Covariances of zooms (from 51 brains).
       0.0009 0.0031 0.0014
       0.0013 0.0014 0.0024];
c44 = diag([0.18 0.11 1.79]*1e-3); % Variance of shears (from 51 brains).
pad = zeros(3);
covar = [c11 pad pad pad
         pad c22 pad pad
         pad pad c33 pad
         pad pad pad c44];
icovar = inv(covar);

ornt = sptl_Ornt;
ornt(7:9) = ornt(7:9).*[1.1 1.05 1.17];

nobayes = 0;

if strcmp(mode,'register1')
	% This is for use in Multimodal coregistration.
	% Each F is mapped to one G (with scaling), but the
	% rigid body components differ between the two sets
	% of registrations.
	%-----------------------------------------------------------------------
	pdesc   = [1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 0
		   0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1]';
	gorder  = [1 2];
	free    = [ones(1,18) ones(1, 2)]';
	mean0   = [ornt([1:6 1:6 7:12]) 1 1]';
	icovar0 = zeros(length(mean0));
	for i=1:size(pdesc,2),
		tmp = find(pdesc(:,i));
		tmp = tmp(1:12);
		icovar0(tmp,tmp) = icovar;
	end;
	ifun = 'spm_matrix([0 0 0  0 0 0  P(7:12)])*spm_matrix(P(1:6))';

elseif strcmp(mode,'rigid1')
	% Rigid body registration.
	% Each F is mapped to one G, without scaling.
	%-----------------------------------------------------------------------
	np     = prod(size(VG));
	if np ~= prod(size(VF))
		error('There should be the same number of object and template images');
	end

	pdesc   = [ones(12,np); eye(np)];
	gorder  = 1:np;
	free    = [ones(1,6) zeros(1,6) zeros(1, np)]';
	mean0   = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
	icovar0 = zeros(length(mean0));
	nobayes = 1;
	ifun = 'spm_matrix(P(1:12))';

elseif strcmp(mode,'rigid2')
	% Rigid body registration.
	% Each F is mapped to one G, with scaling.
	%-----------------------------------------------------------------------
	np     = prod(size(VG));
	if np ~= prod(size(VF))
		error('There should be the same number of object and template images');
	end

	pdesc   = [ones(12,np); eye(np)];
	gorder  = 1:np;
	free    = [ones(1,6) zeros(1,6) ones(1, np)]';
	mean0   = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
	icovar0 = zeros(length(mean0));
	nobayes = 1;
	ifun = 'spm_matrix(P(1:12))';

elseif strcmp(mode,'rigid3')
	% Rigid body registration.
	% Each F is mapped to a linear combination of Gs.
	%-----------------------------------------------------------------------
	np    = prod(size(VG));
	if prod(size(VF)) ~= 1
		error('There should be one object image');
	end

	pdesc   = ones(12+np,1);
	gorder  = ones(1,np);
	free    = [ones(1,6) zeros(1,6) ones(1, np)]';
	mean0   = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
	icovar0 = zeros(length(mean0));
	nobayes = 1;
	ifun = 'spm_matrix(P(1:12))';

elseif strcmp(mode,'2d1')
	% Rigid body registration.
	% Each F is mapped to a linear combination of Gs.
	%-----------------------------------------------------------------------
	np    = prod(size(VG));
	if prod(size(VF)) ~= 1
		error('There should be one object image');
	end
	pdesc   = ones(12+np,1);
	gorder  = ones(1,np);
	free    = [[1 1 0 0 0 1] zeros(1,6) ones(1, np)]';
	mean0   = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
	icovar0 = zeros(length(mean0));
	nobayes = 1;
	ifun = 'spm_matrix(P(1:12))';

elseif strcmp(mode,'affine1')
	% Affine normalisation.
	% Each F is mapped to one G, without scaling.
	%-----------------------------------------------------------------------
	np     = prod(size(VG));
	if np ~= prod(size(VF))
		error('There should be the same number of object and template images');
	end

	pdesc   = [ones(12,np); eye(np)];
	gorder  = 1:np;
	free    = [ones(1,12) zeros(1, np)]';
	mean0   = [ornt ones(1,np)]';
	icovar0 = zeros(length(mean0));
	icovar0(1:12,1:12) = icovar(1:12,1:12);
	ifun = 'spm_matrix(P(1:12))';

elseif strcmp(mode,'affine2')
	% Affine normalisation.
	% Each F is mapped to one G, with scaling.
	%-----------------------------------------------------------------------
	np     = prod(size(VG));
	if np ~= prod(size(VF))
		error('There should be the same number of object and template images');
	end

	pdesc   = [ones(12,np); eye(np)];
	gorder  = 1:np;
	free    = [ones(1,12) ones(1, np)]';
	mean0   = [ornt ones(1,np)]';
	icovar0 = zeros(length(mean0));
	icovar0(1:12,1:12) = icovar(1:12,1:12);
	ifun = 'spm_matrix(P(1:12))';

elseif strcmp(mode,'affine3')
	% Affine normalisation.
	% Each F is mapped to a linear combination of Gs.
	%-----------------------------------------------------------------------
	np    = prod(size(VG));
	if prod(size(VF)) ~= 1
		error('There should be one object image');
	end

	pdesc   = ones(12+np,1);
	gorder  = ones(1,np);
	free    = [ones(1,12) ones(1, np)]';
	mean0   = [ornt ones(1,np)]';
	icovar0 = zeros(length(mean0));
	icovar0(1:12,1:12) = icovar(1:12,1:12);
	ifun = 'spm_matrix(P(1:12))';
else
	error('I don''t understand');
end

if nargin < 6,
	params = mean0;
else,
	% Allow pass of empty params
	if isempty(params),
		params = mean0;
	end;
end;

if nargin<7, VW  = []; end;
if nargin<8, VW2 = []; end;

% Do the optimisation
%-----------------------------------------------------------------------
if nobayes == 1
	[params] = spm_affsub2(ifun,VG,VF,VW,VW2, Hold,samp,params,free,pdesc,gorder);
else
	[params] = spm_affsub2(ifun,VG,VF,VW,VW2, Hold,samp,params,free,pdesc,gorder,mean0,icovar0);
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
function P = spm_affsub2(ifun,VG,VF,VW,VW2, Hold,samp, P,free,pdesc,gorder,mean0,icovar0)
% Another subroutine involved in affine transformations.
% FORMAT P = spm_affsub2(ifun,VG,VF,VW,VW2,Hold,samp,P,free,pdesc,gorder)
%
% ifun      - Function generating affine transformation matrix from 12
%             parameters.
% VG        - Vector of template volumes.
% VF        - Object volume.
% VW        - Memory mapped weighting volume(s) for template image(s)
% VW2       - Memory mapped weighting volume(s) for object image(s)
% Hold      - Interpolation method.
% samp      - Frequency (in mm) of sampling.
% P         - Old parameter estimates.
% free      - Ones and zeros indicating which parameters to fit.
% pdesc     - Description of parameters.
% gorder    - Order in which the template images are used.
% mean0     - The mean of the a-priori probability distribution.
% icovar0   - Inverse of covariance matrix describing the prob. dist.
%
% (Returns) P - New parameter estimates.
%__________________________________________________________________________
%
% Sorry, but the clearest description of what this subroutine does can only
% be obtained by reading the Matlab code.
%__________________________________________________________________________


% Minimal amount of input checking.
%-----------------------------------------------------------------------
if nargin ~= 13 & nargin ~= 11, error('Incorrect usage.'); end;
tmp = sum(pdesc ~= 0);
if size(tmp,2) ~= size(VF,1),
	error(['Incompatible number of object images']);
end;
for i=1:length(tmp),
	if tmp(i) ~= 12+sum(gorder == i),
		error(['Problem with column ' num2str(i) ' of pdesc']);
	end;

end;
if any(gorder > size(VG,1)), error(['Problem with gorder']); end;
if ~all(size(free) == size(P)) | ~(size(pdesc,1) == size(P,1)) | size(P,2) ~= 1,
	error('Problem with vector sizes');
end;

if nargin == 13,
	useW = 1;
	if any(size(mean0) ~= size(P))
		error('A-priori means are wrong size');
	end
	if any(size(icovar0) ~= length(P))
		error('A-priori inv-covariance is wrong size');
	end
else
	useW    = 0;
	mean0   = P;
	icovar0 = diag(eps*ones(prod(size(mean0)),1));
end

iter       = 1;
countdown  = 0;
logdet     = 0;
bestlogdet = 0;

qq    = find(free);
IC0   = icovar0(qq,qq);
P0    = mean0(qq);

W = ones(size(pdesc,2),3)*Inf;
while iter <= 128 & countdown < 4,

	% generate alpha and beta
	%-----------------------------------------------------------------------
	alpha = zeros(length(P));
	beta  = zeros(length(P),1);

	for im = 1:size(pdesc,2),
		pp  = find(pdesc(:,im));
		vf  = VF(im);
		vg  = VG(find(gorder == im));
		if ~isempty(VW ), vw   = VW(im ); else, vw   = []; end;
		if ~isempty(VW2), vw2  = VW2(im); else, vw2  = []; end;

		if useW,
			[alpha_t, beta_t, chi2_t, W(im,:)] = ...
				spm_affsub1(ifun,vg, vf, vw, vw2, Hold,samp,P(pp),W(im,:));
		else,
			[alpha_t, beta_t, chi2_t] = ...
				spm_affsub1(ifun,vg, vf, vw, vw2, Hold,samp,P(pp));
		end;

		beta(pp)     = beta(pp)     + beta_t;
		alpha(pp,pp) = alpha(pp,pp) + alpha_t;
	end;

	% Remove the `fixed' elements
	%----------------------------------------------------------------------
	alpha = alpha(qq,qq);
	beta  = beta(qq);

	% This should give a good indication of the tightness of the fit
	% - providing that the images are all independant.
	% However, future work may involve optimizing Wilk's Lambda for
	% multivariate image registration (i.e., find the registration
	% parameters that maximise the dependance of a linear combination
	% of one set of images upon another set).
	%----------------------------------------------------------------------
	logdet = sum(log(eps+svd(alpha+IC0)));
	spm_chi2_plot('Set', logdet);

	% Check stopping criteria. If satisfied then just do another few more
	% iterations before stopping.
	%-----------------------------------------------------------------------
	if 2*(logdet-bestlogdet)/(logdet+bestlogdet) < 0.0002, countdown = countdown + 1;
	else, countdown = 0; end;

	% If the likelihood is better than the previous best, then save the
	% parameters from the previous iteration.
	%-----------------------------------------------------------------------
	if logdet > bestlogdet & iter > 1, bestlogdet = logdet; end;

	% Update parameter estimates
	%----------------------------------------------------------------------
	P(qq) = pinv(alpha + IC0) * (alpha*P(qq) - beta + IC0*P0);

	iter = iter + 1;
end;
return;
%__________________________________________________________________________

%__________________________________________________________________________
function [alpha, beta, chi2, W] = spm_affsub1(ifun,VG,VF,VW,VW2,Hold,samp,P,minW)
% Generate A'*A and A'*b and \Chi^2 for affine image registration.
% FORMAT [alpha, beta, chi2, W] = spm_affsub1(VG,VF,VW,VW2,Hold,samp,P,minW)
% ifun      - Function generating affine transformation matrix from 12
%             parameters.
% VG        - Vector of template volume(s) (see spm_vol).
% VF        - Object volume.
% VW        - Weighting volume (for template).
% VW2       - Weighting volume (for object).
% Hold      - Interpolation method.
% samp      - frequency (in mm) of sampling.
% P         - Current parameter estimates.
% minW      - previous minimum smoothness estimate.
% 
% alpha     - A'*A
% beta      - A'*b
% chi2      - Residual sum of squares.
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

fun = inline(ifun,'P');

% Flag for masking of template or object
if ~isempty(VW) | ~isempty(VW2)
	wF = 1;
else
	wF = 0;
end

% Sample about every samp mm
%-----------------------------------------------------------------------
vx = sqrt(sum(VG(1).mat(1:3,1:3).^2));
skipx = max([samp/vx(1) 1]);
skipy = max([samp/vx(2) 1]);
skipz = max([samp/vx(3) 1]);


% Convert parameters to affine transformation matrix
%-----------------------------------------------------------------------
Mat = inv(VG(1).mat\fun(P')*VF(1).mat);

% rate of change of matrix elements with respect to parameters
%-----------------------------------------------------------------------
dMdP = zeros(12+length(VG),(12+length(VG)));
tmp  = Mat(1:3,1:4)';
t0   = [tmp(:); zeros(length(VG),1)];
for pp = 1:12;
	tP         = P;
	tP(pp)     = tP(pp)+0.001;
	tmp        = inv(VG(1).mat\fun(tP')*VF(1).mat);
	tmp        = tmp(1:3,1:4)';
	dMdP(:,pp) = ([tmp(:); zeros(length(VG),1)]-t0)/0.001;
end
dMdP(:,(1:length(VG))+12) = [zeros(12,length(VG)); eye(length(VG))];

% Initialise variables
%-----------------------------------------------------------------------
alpha = zeros(12+length(VG),12+length(VG));
beta  = zeros(12+length(VG),1);
chi2  = 0;
dch2  = [0 0 0];
n     = 0;

for p=1:skipz:VG(1).dim(3),	% loop over planes

	% Coordinates of templates
	%-----------------------------------------------------------------------
	[Y,X] = meshgrid(1:skipy:VG(1).dim(2), 1:skipx:VG(1).dim(1));
	X=X(:); Y=Y(:);

	% Transformed template coordinates.
	%-----------------------------------------------------------------------
	X1 = Mat(1,1)*X + Mat(1,2)*Y + (Mat(1,3)*p + Mat(1,4));
	Y1 = Mat(2,1)*X + Mat(2,2)*Y + (Mat(2,3)*p + Mat(2,4));
	Z1 = Mat(3,1)*X + Mat(3,2)*Y + (Mat(3,3)*p + Mat(3,4));

	if wF,
		if ~isempty(VW),
			% Sample weighting image
			%---------------------------------------------------------------
			MatW  = inv(VG(1).mat\VW(1).mat);
			XW    = MatW(1,1)*X + MatW(1,2)*Y + (MatW(1,3)*p + MatW(1,4));
			YW    = MatW(2,1)*X + MatW(2,2)*Y + (MatW(2,3)*p + MatW(2,4));
			ZW    = MatW(3,1)*X + MatW(3,2)*Y + (MatW(3,3)*p + MatW(3,4));
			wt    = spm_sample_vol(VW(1), XW, YW, ZW, 1);
			if ~isempty(VW2),
				MatW  = inv(VG(1).mat\fun(P')*VW2(1).mat);
				XW    = MatW(1,1)*X + MatW(1,2)*Y + (MatW(1,3)*p + MatW(1,4));
				YW    = MatW(2,1)*X + MatW(2,2)*Y + (MatW(2,3)*p + MatW(2,4));
				ZW    = MatW(3,1)*X + MatW(3,2)*Y + (MatW(3,3)*p + MatW(3,4));
				wt2   = spm_sample_vol(VW2(1), XW, YW, ZW, 1);
				wt    = ((wt+eps).^(-1) + (wt2+eps).^(-1)).^(-1);
			end;
		else,    %only object weights
			MatW  = inv(VG(1).mat\fun(P')*VW2(1).mat);
			XW    = MatW(1,1)*X + MatW(1,2)*Y + (MatW(1,3)*p + MatW(1,4));
			YW    = MatW(2,1)*X + MatW(2,2)*Y + (MatW(2,3)*p + MatW(2,4));
			ZW    = MatW(3,1)*X + MatW(3,2)*Y + (MatW(3,3)*p + MatW(3,4));
			wt    = spm_sample_vol(VW2(1), XW, YW, ZW, 1);
		end;

		% Only resample from within the volume VF and where the weight > 0.005.
		%-----------------------------------------------------------------------
		t = 4.9e-2;
		mask1 = find((Z1>=1-t) & (Z1<=VF(1).dim(3)+t) ...
			   & (Y1>=1-t) & (Y1<=VF(1).dim(2)+t) ...
			   & (X1>=1-t) & (X1<=VF(1).dim(1)+t) & wt>0.005);
		wt = sqrt(wt(mask1));
	else,
		% Only resample from within the volume VF.
		%-----------------------------------------------------------------------
		t = 4.9e-2;
		mask1 = find((Z1>=1-t) & (Z1<=VF(1).dim(3)+t) ...
			   & (Y1>=1-t) & (Y1<=VF(1).dim(2)+t) ...
			   & (X1>=1-t) & (X1<=VF(1).dim(1)+t));
	end;


	% Don't waste time on an empty plane.
	%-----------------------------------------------------------------------
	if length(mask1>0),

		% Only resample from within the volume VF.
		%-----------------------------------------------------------------------
		if length(mask1) ~= prod(size(X1)),
			X1 = X1(mask1);
			Y1 = Y1(mask1);
			Z1 = Z1(mask1);
			X  = X(mask1);
			Y  = Y(mask1);
		end;
		Z = zeros(size(mask1))+p;

		% Rate of change of residuals w.r.t parameters
		%-----------------------------------------------------------------------
		dResdM = zeros(size(mask1,1),12+length(VG));

		% Sample object image & get local derivatives
		%-----------------------------------------------------------------------
		[F,dxF,dyF,dzF] = spm_sample_vol(VF, X1, Y1, Z1, Hold);

		% Sample referance image(s) and derivatives
		%-----------------------------------------------------------------------
		for i=1:length(VG),
			if nargout>=4,
				% For computing gradients of residuals
				[Gi,dxt,dyt,dzt] = spm_sample_vol(VG(i), X, Y, Z, Hold);
				if i==1,
					res = F   - Gi*P(i+12); % Residuals
					dxG = dxt*P(i+12);
					dyG = dyt*P(i+12);
					dzG = dzt*P(i+12);
				else,
					res = res - Gi*P(i+12);
					dxG = dxG+dxt*P(i+12);
					dyG = dyG+dyt*P(i+12);
					dzG = dzG+dzt*P(i+12);
				end;
			else,
				Gi = spm_sample_vol(VG(i), X, Y, Z, Hold);
				if i==1, res = F   - Gi*P(i+12); % Residuals
				else,    res = res - Gi*P(i+12); end;
			end;
			if wF, Gi = Gi.*wt; end;
			dResdM(:,12+i) = -Gi;
		end;

		if wF,
			dxF  = dxF.*wt;
			dyF  = dyF.*wt;
			dzF  = dzF.*wt;
			if nargout>=4,
				dxG  = dxG.*wt;
				dyG  = dyG.*wt;
				dzG  = dzG.*wt;
			end;
			res = res.*wt;
		end;

		% Generate Design Matrix from rate of change of residuals wrt matrix
		% elements.
		%-----------------------------------------------------------------------
		dResdM(:,1:12) = [	X.*dxF Y.*dxF p*dxF dxF ...
	       	            		X.*dyF Y.*dyF p*dyF dyF ...
	       	            		X.*dzF Y.*dzF p*dzF dzF ];

		% alpha = alpha + A'*A and beta = beta + A'*b
		%-----------------------------------------------------------------------
		alpha = alpha + spm_atranspa(dResdM);
		beta  = beta  + dResdM'*res;
		clear dResdM

		% Assorted variables which are used later.
		%-----------------------------------------------------------------------
		chi2  = chi2 + res'*res;		% Sum of squares of residuals
		if wF, n = n + sum(wt.*wt);
		else,  n = n + prod(size(F)); end;

		if nargout>=4,
			% Spatial derivatives of residuals derived from
			% (derivatives of F rotated to space of G) - (derivatives of G).
			%-----------------------------------------------------------------------
			tmp = Mat(1:3,1:3)';
			dF  = [dxF dyF dzF]*tmp';
			dxF  = dF(:,1) - dxG;
			dyF  = dF(:,2) - dyG;
			dzF  = dF(:,3) - dzG;
			dch2  = dch2 + [dxF'*dxF dyF'*dyF dzF'*dzF];	% S.o.sq of derivs of residls
		end;
	end;
end;

if n<=32,
	str = {...
		'There is not enough overlap in the',...
		'    images to obtain a solution.',...
		' ',...
		'  Please check that your header information is OK.'};

	spm('alert*',str,mfilename,sqrt(-1));

	error('There is not enough overlap of the images to obtain a solution');
end;

% Smoothness estimate from residuals to determine number of
% independant observations.
%-----------------------------------------------------------------------
if nargout>=4,
	W      = (2*dch2/chi2).^(-.5).*vx;
	msk    = find(~finite(W));
	W      = min([W;minW]);
	W(msk) = 1;
	skips  = [skipx skipy skipz].*vx;		% sample distances (mm)
	skips(msk)  = 1;
	smo    = prod(min(skips./(W*sqrt(2*pi)),[1 1 1])); % fmri revisited
else,
	smo = 1;
end;

df     = (n - size(beta,1))*smo;

% Compute alpha and beta
%-----------------------------------------------------------------------
chi2  = chi2/df;
alpha = dMdP'*alpha*dMdP/chi2;
beta  = dMdP'*beta      /chi2;

return;

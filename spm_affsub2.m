function P = spm_affsub2(VG,VF, MG,MF, Hold,samp, P,mean0,covar,pdesc,gorder)
% Another subroutine involved in affine transformations.
% FORMAT nP = spm_affsub2(VG,VF,MG,MF,Hold,samp,oP,covar,mean0,pdesc,gorder)
%
% VG        - Vector of memory mapped template image(s).
% VF        - Memory mapped object image.
% MG        - Space of the template image(s).
% MF        - Space of the object image.
% Hold      - Interpolation method.
% samp      - Frequency (in mm) of sampling.
% oP        - Old parameter estimates.
% mean0     - A-priori mean value of parameters.
% covar     - A-priori covariance matrix of parameters.
% pdesc     - Description of parameters.
% gorder    - Order in which the template images are used.
%
% nP        - New parameter estimates.
%__________________________________________________________________________
%
% Sorry, but the clearest description of what this subroutine does can only
% be obtained by reading the Matlab code.
%__________________________________________________________________________
% %W% John Ashburner FIL %E%


% Minimal amount of input checking.
%-----------------------------------------------------------------------
if nargin ~= 11
	error('Incorrect usage.');
end
if size(VF,2) ~= size(MF,2) | size(VG,2) ~= size(MG,2)
	error('Incompatible number of position matrixes');
end
tmp = sum(pdesc ~= 0);
if size(tmp,2) ~= size(VF,2)
	error(['Incompatible number of object images']);
end
for i=1:length(tmp)
	if tmp(i) ~= 12+sum(gorder == i)
		error(['Problem with column ' num2str(i) ' of pdesc']);
	end

end
if any(gorder > size(VG,2))
	error(['Problem with gorder']);
end
if ~all([size(covar) size(mean0,1)] == size(pdesc,1)) | size(mean0,2) ~= 1
	error('Problem with covar or mean0');
end


% Updating scheme is a weighted average of mean0 and (P + alpha\beta):
% P = (alpha + inv(covar))\(alpha\(P + alpha\beta) + covar\mean0)
%-----------------------------------------------------------------------
alpha0    = pinv(covar);
beta0     = alpha0*mean0;

pchi2     = 9e99;
iter      = 1;
countdown = 0;
bestP     = P;
bestchi2  = 9e99;

while iter <= 64 & countdown < 3

	ochi2 = pchi2;
	pchi2 = 1;

	% generate alpha and beta
	%-----------------------------------------------------------------------
	alpha = alpha0;
	beta  = beta0;
	for im = 1:size(pdesc,2)

		pp  = find(pdesc(:,im));
		ppp = find(pdesc(:,im)*pdesc(:,im)');

		vf  = VF(:,im);
		vg  = VG(:,find(gorder == im));
		mf  = reshape(MF(:,im),4,4);

		% Note: only the matrix of the first template image is used
		%       when registering one image to many.
		%-----------------------------------------------------------------------
		tmp = find(gorder == im);
		tmp = tmp(1);
		mg  = reshape(MG(:,tmp),4,4);

		[alpha_t, beta_t, chi2_t] = spm_affsub1(vg, vf, mg, mf, Hold,samp,P(pp));

		beta(pp)   = beta(pp)   + beta_t(:) /chi2_t;
		alpha(ppp) = alpha(ppp) + alpha_t(:)/chi2_t;
		pchi2      = pchi2 * chi2_t;
	end

	% If \chi^2 is better than the previous best, then save the parameters
	% from the previous iteration.
	%-----------------------------------------------------------------------
	if (pchi2 < bestchi2)
		bestchi2  = pchi2;
		bestP     = P;
	end

	% Update parameter estimates
	%-----------------------------------------------------------------------
	P = pinv(alpha)*(alpha*P + beta);

	% Check stopping criteria. If satisfied then just do another few more
	% iterations before stopping.
	%-----------------------------------------------------------------------
	if (2*(ochi2-pchi2)/(ochi2 + pchi2)) < 0.01
		countdown = countdown + 1;
	end

	iter = iter + 1;
end

P = bestP;

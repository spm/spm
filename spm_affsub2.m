function P = spm_affsub2(VG,VF, MG,MF, Hold,samp, P,free,pdesc,gorder)
% Another subroutine involved in affine transformations.
% FORMAT nP = spm_affsub2(VG,VF,MG,MF,Hold,samp,oP,free,pdesc,gorder)
%
% VG        - Vector of memory mapped template image(s).
% VF        - Memory mapped object image.
% MG        - Space of the template image(s).
% MF        - Space of the object image.
% Hold      - Interpolation method.
% samp      - Frequency (in mm) of sampling.
% oP        - Old parameter estimates.
% free      - Ones and zeros indicating which parameters to fit.
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
if nargin ~= 10
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
if ~all(size(free) == size(P)) | ~(size(pdesc,1) == size(P,1)) | size(P,2) ~= 1
	error('Problem with vector sizes');
end

pchi2     = 9e99;
iter      = 1;
countdown = 0;
bestP     = P;
bestchi2  = 9e99;

qq  = find(free);
qqq = find(free*free');
nf  = sum(free ~= 0);

while iter <= 64 & countdown < 3

	ochi2 = pchi2;
	pchi2 = 1;

	% generate alpha and beta
	%-----------------------------------------------------------------------
	alpha = zeros(length(P));
	beta  = zeros(length(P),1);

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
	P(qq) = P(qq) + pinv(reshape(alpha(qqq), nf, nf))*beta(qq);

	% Check stopping criteria. If satisfied then just do another few more
	% iterations before stopping.
	%-----------------------------------------------------------------------
	if (2*(ochi2-pchi2)/(ochi2 + pchi2)) < 0.01
		countdown = countdown + 1;
	end

	iter = iter + 1;
end

P = bestP;

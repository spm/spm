function [P] = spm_affsub2(VG,VF,VW, Hold,samp, P,free,pdesc,gorder,mean0,icovar0)
% Another subroutine involved in affine transformations.
% FORMAT [nP] = spm_affsub2(VG,VF,VW,Hold,samp,oP,free,pdesc,gorder)
%
% VG        - Vector of memory mapped template image(s).
% VF        - Memory mapped object image.
% Hold      - Interpolation method.
% samp      - Frequency (in mm) of sampling.
% oP        - Old parameter estimates.
% free      - Ones and zeros indicating which parameters to fit.
% pdesc     - Description of parameters.
% gorder    - Order in which the template images are used.
% mean0     - The mean of the a-priori probability distribution.
% icovar0   - Inverse of covariance matrix describing the prob. dist.
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
if nargin ~= 11 & nargin ~= 9
	error('Incorrect usage.');
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

if nargin == 12
	useW=1;
	if any(size(mean0) ~= size(P))
		error('A-priori means are wrong size');
	end
	if any(size(icovar0) ~= length(P))
		error('A-priori inv-covariance is wrong size');
	end
else
	useW=0;
	mean0 = P;
	icovar0 = diag(eps*ones(prod(size(mean0)),1));
end

iter       = 1;
countdown  = 0;
bestP      = P;
logdet     = 0;
bestlogdet = 0;

qq    = find(free);
qqq   = find(free*free');
nf    = sum(free ~= 0);
IC0   = reshape(icovar0(qqq), nf, nf);
P0    = mean0(qq);

W = ones(size(pdesc,2),3)*Inf;
while iter <= 128 & countdown < 4


	% generate alpha and beta
	%-----------------------------------------------------------------------
	alpha = zeros(length(P));
	beta  = zeros(length(P),1);

	for im = 1:size(pdesc,2)

		pp  = find(pdesc(:,im));
		ppp = find(pdesc(:,im)*pdesc(:,im)');

		vf  = VF(im);
		vg  = VG(find(gorder == im));
		if ~isempty(VW)
			vw  = VW(im);
		else
			vw  = [];
		end

		% Note: only the matrix of the first template image is used
		%       when registering one image to many.
		%-----------------------------------------------------------------------
		tmp = find(gorder == im);
		tmp = tmp(1);

		%if iter>1 flg=1; else flg=0; end
		flg=0;

		if useW,
			[alpha_t, beta_t, chi2_t, W(im,:)] = ...
				spm_affsub1(vg, vf, vw, Hold,samp,P(pp),flg,W(im,:));
		else
			[alpha_t, beta_t, chi2_t] = ...
				spm_affsub1(vg, vf, vw, Hold,samp,P(pp),flg);
		end

		beta(pp)   = beta(pp)   + beta_t(:);
		alpha(ppp) = alpha(ppp) + alpha_t(:);
	end

	% Remove the `fixed' elements
	%----------------------------------------------------------------------
	alpha = reshape(alpha(qqq), nf, nf);
	beta  = beta(qq);

	% This should give a good indication of the tightness of the fit
	%----------------------------------------------------------------------
	logdet = sum(log(eps+svd(alpha+IC0)));

	spm_chi2_plot('Set', logdet);
	%fprintf('iteration: %d\tlog(det): %g\n', iter, logdet);

	% Check stopping criteria. If satisfied then just do another few more
	% iterations before stopping.
	%-----------------------------------------------------------------------
	if (2*(logdet-bestlogdet)/(logdet+bestlogdet) < 0.0002) countdown = countdown + 1;
	else countdown = 0; end;

	% If the likelyhood is better than the previous best, then save the
	% parameters from the previous iteration.
	%-----------------------------------------------------------------------
	if (logdet > bestlogdet & iter > 1)
		bestlogdet = logdet;
		bestP  = P;
	end;

	% Update parameter estimates
	%----------------------------------------------------------------------
	P(qq) = pinv(alpha + IC0) * (alpha*P(qq) - beta + IC0*P0);

	iter = iter + 1;
end

%P = bestP;

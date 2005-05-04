function [w,ind,b0] = spm_rvc(K,t,bs,alpha,th)
% Optimisation for Relevance Vector Classification
% USAGE [w,ind,b0] = spm_rvc(K,t,bs)
% w    - non-zero weights 
% ind  - index for non-zero elements
% K    - MxN matrix derived from kernel function of vector pairs
% t    - M vector showing group memberships (1 or 0)
% bs   - include bias term
%
% The REML solution tends towards infinite weights for some the
% regularisation terms (i.e. 1/alpha(i) approaches 0).
% The appropriate columns are removed from the model when
% this happens.
%
% see: http://research.microsoft.com/mlp/RVM/relevance.htm
%
% Refs:
% The Relevance Vector Machine.
% In S. A. Solla, T. K. Leen, and K.-R. Müller (Eds.),
% Advances in Neural Information Processing Systems 12,
% pp.  652-658. Cambridge, Mass: MIT Press.
%
% Michael E. Tipping
% Sparse Bayesian Learning and the Relevance Vector Machine
% Journal of Machine Learning Research 1 (2001) 211-244
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_rvc.m 112 2005-05-04 18:20:52Z john $


if nargin<3, bs = 0; end;
if nargin<5,
	s = max(svd(K))^2;
	if nargin<4, alpha = s*max(size(K))*1e-10; end;
	if nargin<5, th    = s*1e8;               end;
end;

if bs,
	K = [K ones(size(K,1),1)];
end;
[N,M]	= size(K);
t       = logical(t);
w	= zeros(M,1);
alpha   = ones(M,1)*alpha;
for i=1:1000
	nz        = alpha<th;
	if sum(nz)==0,
		warning('Converged to an empty solution: using previous weights.');
		nz  = alpha_old<th;
		ind = find(nz);
		w   = w(ind);
		break;
	end;
	alpha_old = alpha;

	% E-step
	[w(nz),C] = map_soln(K(:,nz),t,alpha(nz),w(nz));

	% M-step
	alpha(nz) = (1 - alpha(nz).*diag(C)) ./ w(nz).^2;

	% % The alternative scheme - but I don't understand its motivation
	% gamma     = (1 - alpha(nz).*diag(C));
	% alpha(nz) = gamma./(w(nz).^2./gamma - diag(C));
	% alpha(alpha<=0) = Inf;

	% Convergence
	if max(abs(log((alpha(nz)+eps)./(alpha_old(nz)+eps)))) < 1e-7,
		break;
	end;
end;
b0  = 0.0;
ind = find(nz);
w   = w(ind);
if bs & ind(end)>N,
	b0  = w(end);
	w   = w(1:(end-1));
	ind = ind(1:(end-1));
end;
return;
%__________________________________________________________________________
 
%__________________________________________________________________________
function [w,C] = map_soln(K,t,alpha,w)
% Levenberg-Marquardt optimisation of MAP solution of w
% Also returns the formal covariance matrix of the fit on the assumption
% of normally distributed errors (inverse of Hessian), for Laplace
% approximation. The main reason for the regularisation is that
% the Hessian can be badly conditioned.

[M,N]   = size(K);
[err_old,y] = errors(K,t,alpha,w);
[g,H]   = grad_hess(K,t,y,alpha,w);
lambda  = max(eig(H))*eps*size(H,1);
I       = eye(size(H));
for it=1:100,
	w_new       = w + (H + I*lambda)\g;
	[err_new,y] = errors(K,t,alpha,w_new);
	if err_new <= err_old,
		err_old = err_new;
		w       = w_new;
		[g,H]   = grad_hess(K,t,y,alpha,w);
		if it>1 & norm(g)<N*1e-7,
			break;
		end
		% With regular LM, the following would be
		% included.
		% lambda = lambda/10;
	else,
		lambda = lambda*10;
	end;
end;
C = inv(H+I*lambda);
return;
%__________________________________________________________________________
 
%__________________________________________________________________________
function [err,y] = errors(K,t,alpha,w)
% Return errors and sigmoid
[M,N] = size(K);
y     = 1./(1+exp(-(K*w)));
err   = -(sum(log(y(t)+eps))+sum(log(1-y(~t)+eps)))/M + (alpha'*(w.^2))/(2*M);
return;
%__________________________________________________________________________
 
%__________________________________________________________________________
function [g,H] = grad_hess(K,t,y,alpha,w)
% Gradient and Hessian of above function
[M,N] = size(K);
g     = K'*(t-y) - alpha.*w;
H     = (K.*repmat(y.*(1-y),1,N))'*K + diag(alpha);
return;
%__________________________________________________________________________

function [Ep,Cp,Ce] = spm_nlsi_GN(M,U,Y,Ce)
% Bayesian Parameter estimation using the Gauss-Newton method/EM algorithm
% FORMAT [Ep,Cp,Ce] = spm_nlsi_GN(M,U,Y,Ce)
%
% M   - model structure
% U   - input  structure
% Y   - output structure
% Ce  - constraints on Cov{e}
%
% Model Parameter estimates
%---------------------------------------------------------------------------
% Ep  - (p x 1)        	conditional expectation  E{P|y}
% Cp  - (p x p)        	conditional covariance   Cov{P|y}
% Ce  - (v x v)        	ML estimate of           Cov{e}
%
%___________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a 
% nonlinear observation model under Gaussian assumptions
%
%  		 	dx/dt = f(x,u,P)
%   		 	y(t)  = l(x(t),P) + e
%
% with priors on the free parameters P specified in terms of expectation pE
% and covariance pC. The estimation uses a Gauss-Newton method with MAP
% point estimators at each iteration.  The covariance of e is a maximum
% likelihood estimator based on the residuals (c.f. the EM algorithm).
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% repsonse variable and confounds if applicable
%---------------------------------------------------------------------------
y      = Y.y;
[v u]  = size(y);
if isfield(Y,'X0')
	X0  = kron(speye(u,u),Y.X0);
else
	X0  = [];
end

% SVD of prior covariances
%---------------------------------------------------------------------------
V      = spm_svd(M.pC,1e-16);
pE     = V'*M.pE(:);
pC     = V'*M.pC*V;

% treat confounds as fixed effects
%---------------------------------------------------------------------------
u      = size(X0,2);
E0     = sparse(u,1);
C0     = speye(u,u)*1e8;
pC     = blkdiag(pC,C0);
P{1}.C = Ce;
P{2}.C = pC;

% Gauss-Newton search 
%===========================================================================
p      = M.pE(:);
dv     = 1e-4;
for  j = 1:32

	% f(p) - for new expansion point (p)
	%-------------------------------------------------------------------
	fx     = spm_nlsi_int(p,M,U,v);

	% compute partial derivatives df(x)/dx {Jacobian}
	%-------------------------------------------------------------------
	for  i = 1:size(V,2)
		pi     = p + V(:,i)*dv;
		dfx    = spm_nlsi_int(pi,M,U,v) - fx;
		J(:,i) = dfx(:)/dv;
	end

	% Bayesian estimator of new expansion point
	%-------------------------------------------------------------------
	P{1}.X = [J X0];
	P{2}.X = [(pE - V'*p); E0];
	[C P]  = spm_PEB(y(:) - fx(:),P);
	dp     = V*C{2}.E(1:size(V,2));
	p      = p + dp;

	% convergence
	%-------------------------------------------------------------------
	w      = dp'*dp;
	fprintf('%-30s: %i %20s%e\n','GNS Iteration',j,'...',full(w));
	if w < 1e-6, break, end
	
end

% rotate back into native parameter space
%---------------------------------------------------------------------------
Ep     = p;
Cp     = V*C{2}.C(1:size(V,2),1:size(V,2))*V';
Ce     = C{1}.M;

function [ p,Cp,Ce] = spm_nlsi_GN(M,U,Y)
% Bayesian Parameter estimation using the Gauss-Newton method/EM algorithm
% FORMAT [Ep,Cp,Ce] = spm_nlsi_GN(M,U,Y)
%
% Dynamical MIMO models
%___________________________________________________________________________
% M.fx - f(x,u,P)
% M.lx - l(x,P)
%	x - states variables
%	u - inputs or causes
%	P - free parameters
% M.pE - prior expectation  - E{P}
% M.pC - prior covariance   - Cov{P}
%
%
% U.u  - inputs
% U.dt - sampline interval
%
% Y.y  - outputs
% Y.X0 - Confounds or null space
% Y.dt - sampling interval for outputs
% Y.Ce - error covariance contraints
% Y.pC - prior covariance of timing error Cov{q}, E{q} = 0 [default = 0]
%
% Static nonlinear models
%___________________________________________________________________________
% M.fp - f(P,u)
%	u - inputs or causes (c.f. design matrix)
%	P - free parameters
% M.pE - prior expectation  - E{P}
% M.pC - prior covariance   - Cov{P}
% U    - u 
% Y    - output structure
%
% Parameter estimates
%---------------------------------------------------------------------------
% Ep  - (p x 1)        	conditional expectation  E{P|y}
% Cp  - (p x p)        	conditional covariance   Cov{P|y}
% Ce  - (v x v)        	[Re]ML estimate of       Cov{e}
%
%___________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a 
% dynamical MIMO input-state-ouput model under Gaussian assumptions
%
%  		 	dx/dt = f(x,u,P)
%   		 	y(t)  = l(x(t),P) + e
%
% if inputs are specified, otherwise static a nonlinear observation 
% model with fixed input or causes u
%
%                       y(t)  = f(P,u) + e
%
% Priors on the free parameters P specified in terms of expectation pE
% and covariance pC. The estimation uses a Gauss-Newton method with MAP
% point estimators at each iteration.  The covariance of e is a maximum
% likelihood estimator based on the residuals (c.f. the ReML algorithm).
% This corresponds to a Gauss-Newton ascent on the conditional probabilty
% p{P|y}
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% Covariance constraints on observation error (if not specified)
%---------------------------------------------------------------------------
y      = Y.y;
[v m]  = size(y);
if isfield(Y,'Ce')
	Ce = Y.Ce;
else
	Ce = spm_Ce(v*ones(1,m));
end

% Covariance constraints on priors
%---------------------------------------------------------------------------
pE     = M.pE;
pC     = M.pC;
if iscell(pC)
	pC = pC{:};
	EB = 1;
else
	EB = 0;
end

% Covariance priors on timiing error (if not specified)
%---------------------------------------------------------------------------
if isfield(Y,'pC')
	qE = sparse(m,1);
	qC = Y.pC;
else
	qE = sparse(m,1);
	qC = sparse(m,m);
end

% confounds (if specified)
%---------------------------------------------------------------------------
if isfield(Y,'X0')
	Ju = kron(speye(m,m),Y.X0);
else
	Ju = [];
end

% treat confounds Ju as fixed effects
%---------------------------------------------------------------------------
u      = size(Ju,2);
uE     = sparse(u,1);
uC     = speye(u,u)*1e+8;

% SVD of prior covariances
%---------------------------------------------------------------------------
Vp     = spm_svd(pC,1e-16);
Vq     = spm_svd(qC,1e-16);
ip     = [1:size(Vp,2)];
iq     = [1:size(Vq,2)] + ip(end);


% If EB: i.e. prior covariance hyperparameter estimation
%---------------------------------------------------------------------------
pC     = Vp'*pC*Vp;
qC     = Vq'*qC*Vq;
if EB
	P{1}.C = Ce;
	P{2}.C = {blkdiag(0*pC, qC, uC), blkdiag(pC, 0*qC, 0*uC)};
	P{3}.C = 1e-8;
	P{3}.X = 1;
else
	P{1}.C = Ce;
	P{2}.C = blkdiag(pC, qC, uC);
end

% Gauss-Newton search 
%===========================================================================
p      = pE;
q      = qE;
dv     = 1e-6;
for  j = 1:32

	% y(t) = f(p) - for new expansion point (p) in parameter space 
	%-------------------------------------------------------------------
	if isfield(M,'fp')
		fp      = feval(M.fp,p,U);
	else
		[fp,fq] = spm_int(p,M,U,v);
		fp      = fp + fq*diag(q);
	end

	% compute partial derivatives [Jp] dy(t)/dp*Vp {p = parameters}
	%-------------------------------------------------------------------
	for  i = 1:size(Vp,2)

		pi      = p + Vp(:,i)*dv;

		if isfield(M,'fp')
			dfp       = feval(M.fp,pi,U) - fp;
		else
			[fpi fqi] = spm_int(pi,M,U,v);
			fpi       = fpi + fqi*diag(q);
			dfp       = fpi - fp;
		end

		Jp(:,i) = dfp(:)/dv;
	end

	% compute partial derivatives [Jq] dy(t)/dq*Vq {q = delay}
	%-------------------------------------------------------------------
	Jq    = [];
	if size(Vq,2)
		for i = 1:m
			Jq = blkdiag(Jq,fq(:,i));
		end
		Jq    = Jq*Vq;
	end

	% Bayesian [conditional] estimator of new expansion point E{p|y}
	%-------------------------------------------------------------------
	P{1}.X = [Jp Jq Ju];
	P{2}.X = [Vp'*(pE - p); Vq'*(qE - q); uE];
	[C P]  = spm_PEB(y(:) - fp(:),P);

	% update
	%-------------------------------------------------------------------
	dp     = Vp*C{2}.E(ip);
	dq     = Vq*C{2}.E(iq);
	p      = p + dp;
	if size(Vq,2), q = q + dq; end


	% convergence
	%-------------------------------------------------------------------
	w      = dp'*dp;
	fprintf('%-30s: %i %20s%e\n','GNS Iteration',j,'...',full(w));
	if w < 1e-6, break, end

	% graphics
	%-------------------------------------------------------------------
	if length(dbstack) < 3
		bar(p)
		xlabel('parameter')
		ylabel('conditional expectation')
		title(sprintf('%s: %i','GNS Iteration',j))
		grid on
		drawnow
	end

end

% outputs
%---------------------------------------------------------------------------
Ep     = p;
Cp     = Vp*C{2}.C(ip,ip)*Vp';
Ce     = C{1}.M;

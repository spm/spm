function varargout = spm_nlsi(M,U,Y)
% nonlinear system identification of a MIMO system
% FORMAT [Ep,Cp,Ce,K0,K1,K2,M0,M1,L] = spm_nlsi(M,U,Y)
% FORMAT [K0,K1,K2,M0,M1,L] = spm_nlsi(M)
%
% Model specification
%---------------------------------------------------------------------------
% M.fx  - dx/dt = f(x,u,P) 	{function string or m-file}
% M.lx  - y(t)  = l(x,P)    	{function string or m-file}
%
% M.pE  - (p x 1)		Prior expectation of p model parameters
% M.pC  - (p x p)		Prior covariance for p model parameters
%
% M.x   - (n x 1)		intial state x(0)
% M.m   - m			number of inputs
% M.n   - n			number of states
% M.l   - l			number of outputs
% M.N	-                       kernel depth
% M.dt	-                       kernel resolution {secs}
%
% System inputs
%---------------------------------------------------------------------------
% U.u   - (u x m)		m inputs
% U.dt  - 			sampling interval for inputs
%
% System outputs
%---------------------------------------------------------------------------
% Y.y   - (v x l)         	l outputs
% Y.X0  - (v x c)		Confounds or null space
% Y.dt  - 			sampling interval for outputs
% Y.Ce  -			obervation error covariance constraints
%
% Model Parameter estimates - conditional moments
%---------------------------------------------------------------------------
% Ep    - (p x 1)         	conditional expectation  E{P|y}
% Cp    - (p x p)        	conditional covariance   Cov{P|y}
% Ce    - (v x v)        	ML estimate of           Cov{e}
%
% System identification     - Volterra kernels
%---------------------------------------------------------------------------
% K0    - (l x 1)             = k0(t)         = y(t)
% K1    - (N x l x m)         = k1i(t,s1)     = dy(t)/dui(t - s1)
% K2    - (N x N x l x m x m) = k2ij(t,s1,s2) = d2y(t)/dui(t - s1)duj(t - s2)
%
% System identification     - Bilinear approximation
%---------------------------------------------------------------------------
% M0    - (n x n)     		dq/dx
% M1    - (n x n x m) 		d2q/dxdu
% L     - (n x l)     		dl/dq
%
%___________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a 
% nonlinear MIMO model under Gaussian assumptions
%
%              dx/dt  = f(x,u,P)
%                y(t) = l(x,P) + e                               (1
%
% evaluated at x(0) = x0, using a Bayesian estimation scheme with priors
% on the model parameters P, specified in terms of expectations and 
% covariance. The estimation uses a Gauss-Newton method with MAP point 
% estimators at each iteration.  Both Volterra kernel and state-space 
% representations of the Bilinear approximation are provided,
% The Bilinear approximation to (1), evaluated at x(0) = x and u = 0 is:
%
%		dq/dt = M0*q + u*M1*q
%		 y(t) = L*q
%
%		 q(t) = [1; x(t) - x(0); kron(x(t) - x(0),x(t) - x(0))]
%
% The associated kernels are derived using closed form expressions
%
%---------------------------------------------------------------------------
% If the inputs U and outputs Y are not specified the model is simply
% characterised in terms of its Volterra kernels and Bilinear
% approximation expanding around M.pE
%
% see also
% spm_nlsi_GN:   Bayesian parameter estimation using an EM/Gauss-Newton method
% spm_bi_reduce: Reduction of a fully nonlinear MIMO system to Bilinear form
% spm_kernels:   Returns global Volterra kernels for a MIMO Bilinear system
%
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% Expansion point (in parameter space) for Bilinear and kernel representations
%---------------------------------------------------------------------------
if nargin == 3

	% get constraints on Cov{e} - C1
	%-------------------------------------------------------------------
	if ~isfield(Y,'Ce')
		[v l] = size(Y.y);
		Ce    = spm_Ce(v*ones(1,l));
	else
		Ce    = Y.Ce;
	end

	% Gauss-Newton/Bayesian/EM estimation
	%===================================================================
	[Ep,Cp,Ce] = spm_nlsi_GN(M,U,Y);

	if nargout < 4, varargout = {Ep,Cp,Ce}; return, end

else
	% Use prior expectation to expand around
	%-------------------------------------------------------------------
	Ep         = M.pE;
end


% Bilinear representation
%===========================================================================
[M0,M1,L]  = spm_bi_reduce(M,Ep);


% Volterra kernels
%===========================================================================
[K0,K1,K2] = spm_kernels(M0,M1,L,M.N,M.dt);


% output arguments
%---------------------------------------------------------------------------
if nargin == 3
	varargout = {Ep,Cp,Ce,K0,K1,K2,M0,M1,L};
else
	varargout = {K0,K1,K2,M0,M1,L};
end

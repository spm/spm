function varargout = spm_nlsi(M,U,Y)
% nonlinear system identification of a MIMO system
% FORMAT [Ep,Cp,Ce,K0,K1,K2,M0,M1,L1,L2] = spm_nlsi(M,U,Y)
% FORMAT [K0,K1,K2,M0,M1,L1,L2]          = spm_nlsi(M)
%
% Model specification
%---------------------------------------------------------------------------
% M.f   - dx/dt = f(x,u,P) 	{function string or m-file}
% M.g   - y     = g(x,u,P)    	{function string or m-file}
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
% Ce    - (v x v)        	ReML estimate of           Cov{e}
%
% System identification     - Volterra kernels
%---------------------------------------------------------------------------
% K0    - (l x 1)             = k0(t)         = y(t)
% K1    - (N x l x m)         = k1i(t,s1)     = dy(t)/dui(t - s1)
% K2    - (N x N x l x m x m) = k2ij(t,s1,s2) = d2y(t)/dui(t - s1)duj(t - s2)
%
% System identification     - Bilinear approximation
%---------------------------------------------------------------------------
% M0    - (n x n)     		df/dq
% M1    - {m}(n x n) 		d2f/dqdu
% L1    - (l x n)     		dg/dq
% L2    - {l}(n x n)     	d2g/dqdq
%
%___________________________________________________________________________
%
% Returns the moments of the posterior p.d.f. of the parameters of a 
% nonlinear MIMO model under Gaussian assumptions
%
%              dx/dt  = f(x,u,P)
%                y    = g(x,u,P) + e                               (1
%
% evaluated at x(0) = x0, using a Bayesian estimation scheme with priors
% on the model parameters P, specified in terms of expectations and 
% covariance. The estimation uses a Gauss-Newton method with MAP point 
% estimators at each iteration.  Both Volterra kernel and state-space 
% representations of the Bilinear approximation are provided,
% The Bilinear approximation to (1), evaluated at x(0) = x and u = 0 is:
%
%		dq/dt = M0*q + u(1)*M1{1}*q + u(2)*M1{2}*q + ....
%		 y(i) = L1(i,:)*q + q'*L2{i}*q;
%
% whwew the states are augmented with a constant
%
%		 q(t) = [1; x(t) - x(0)]
%
% The associated kernels are derived using closed form expressions based
% on the bilinear approximation.
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
% SEE NOTES AT THE END OF THIS SCRIPT FOT EXAMPLES
%
%---------------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$


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
[M0,M1,L1,L2] = spm_bi_reduce(M,Ep);


% Volterra kernels
%===========================================================================

% time bins (if not specified)
%---------------------------------------------------------------------------
try 
	dt   = M.dt;
	N    = M.N;
catch
	s    = real(eig(full(M0)));
	s    = max(s(find(s < 0)));
	N    = 32;
	dt   = -4/(s*N);
	M.dt = dt;
	M.N  = N;
end

% get kernels
%---------------------------------------------------------------------------
[K0,K1,K2]   = spm_kernels(M0,M1,L1,L2,N,dt);

% graphics
%===========================================================================
if length(dbstack) < 2

	subplot(2,1,1)
	plot([1:N]*dt,K1(:,:,1))
	xlabel('time')
	ylabel('response')
	title('1st-order kernel')
	grid on

	subplot(2,1,2)
	imagesc([1:N]*dt,[1:N]*dt,K2(:,:,1,1,1))
	xlabel('time')
	ylabel('time')
	title('2nd-order kernel')
	grid on
	axis image
	drawnow
end


% output arguments
%---------------------------------------------------------------------------
if nargin == 3
	varargout = {Ep,Cp,Ce,K0,K1,K2,M0,M1,L1,L2};
else
	varargout = {K0,K1,K2,M0,M1,L1,L2};
end


return

% NOTES ON USE
%===========================================================================
% Consider the system
%
%              dx/dt  = 1./(1 + exp(-x*P)) + u
%                y    = x + e
%
% Specify a model structure:
%---------------------------------------------------------------------------
M.f  = inline('1./(1 + exp(-P*x)) + [u; 0]','x','u','P');
M.g  = inline('x','x','u','P');
M.pE = [-1 .3;.5 -1];			% Prior expectation of parameters
M.pC = speye(4,4);			% Prior covariance for parameters
M.x  = zeros(2,1)			% intial state x(0)
M.m  = 1;				% number of inputs
M.n  = 2;				% number of states
M.l  = 2;				% number of outputs

% characterise M in terms of Volterra kernels and Bilinear matrices:
%---------------------------------------------------------------------------
[K0,K1,K2,M0,M1,L1,L2] = spm_nlsi(M);


% or estimate the model parameters with inputs and outputs
%===========================================================================

% inputs
%---------------------------------------------------------------------------
U.name = 'input'
U.u    = randn(128,1);
U.dt   = 1;

% outputs
%---------------------------------------------------------------------------
Y.name = 'response';
y      = spm_int(M.pE,M,U,64);
Y.y    = y + randn(size(y))/32;
Y.dt   = U.dt*length(U.u)/length(Y.y);

% estimate
%---------------------------------------------------------------------------
[Ep,Cp,Ce,K0,K1,K2,M0,M1,L1,L2] = spm_nlsi(M,U,Y);

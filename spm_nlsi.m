function varargout = spm_nlsi(M,U,Y)
% nonlinear system identification of a MIMO system
% FORMAT [Ep,Cp,Ce,H0,H1,H2,K0,K1,K2,A,B,C,D] = spm_nlsi(M,U,Y)
% FORMAT          [H0,H1,H2,K0,K1,K2,A,B,C,D] = spm_nlsi(M)
%
% Model specification
%---------------------------------------------------------------------------
% M.fx  - dx/dt = f(x,u,P) 	{function string or m-file}
% M.lx  - y(t)  = l(x,P)    	{function string or m-file}
%
% M.pE  - (p x 1)		Prior expectation of p model parameters
% M.pC  - (p x p)		Prior covariance for p model parameters
%
% M.x   - (n x 1)		intial state variables x(0)
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
%
% Model Parameter estimates - conditional moments
%---------------------------------------------------------------------------
% Ep    - (p x 1)         	conditional expectation  E{P|y}
% Cp    - (p x p)        	conditional covariance   Cov{P|y}
% Ce    - (v x v)        	ML estimate of           Cov{e}
%
% System identification     - Volterra kernels - state variables
%---------------------------------------------------------------------------
% H0    - (n x 1)             = h0(t)         = x(t)
% H1    - (N x n x m)         = h1i(t,s1)     = dx(t)/dui(t - s1)
% H2    - (N x N x n x m x m) = h2ij(t,s1,s2) = d2x(t)/dui(t - s1)duj(t - s2)
%
% System identification     - Volterra kernels - ouputs
%---------------------------------------------------------------------------
% K0    - (l x 1)             = k0(t)         = y(t)
% K1    - (N x l x m)         = k1i(t,s1)     = dy(t)/dui(t - s1)
% K2    - (N x N x l x m x m) = k2ij(t,s1,s2) = d2y(t)/dui(t - s1)duj(t - s2)
%
% System identification     - Bilinear approximation
%---------------------------------------------------------------------------
% A     - (n x n)     		df(x(0),0)/dx
% B     - (n x n x m) 		d2f(x(0),0)/dxdu
% C     - (n x m)     		df(x(0),0)/du - d2f(x(0),0)/dxdu*x(0)
% D     - (n x 1)     		f(x(0).0) - df(x(0),0)/dx*x(0)
%
%___________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a 
% nonlinear MIMO model under Gaussian assumptions
%
%                           dx/dt  = f(x,u,P)
%                             y(t) = l(x,P) + e                       (1
%
% evaluated at x(0) = x0, using a Bayesian estimation scheme with priors
% on the model parameters P, specified in terms of expectations and 
% covariance. The estimation uses a Gauss-Newton method with MAP point 
% estimators at each iteration.  Both Volterra kernel and state-space 
% representations of the Bilinear approximation are provided, where
%
%   dx/dt = f(x,u) = A*x + B1*x*u1 + ... Bm*x*um + C1u1 + ... Cmum + D
%    y(t) = l(x(t))
%
% If the inputs U and outputs Y are not specified the model is simply
% characterised in terms of its Volterra kernels and Bilinear
% approximation expanding around M.pE
%
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% dimensions
%---------------------------------------------------------------------------
n          = M.n;				% m inputs
m          = M.m;				% l outputs
l          = M.l;				% n state variables

% Determine expansion point for Bilinear and kernel representations
%---------------------------------------------------------------------------
if nargin == 3

	% get dimensions of inputs and outputs
	%-------------------------------------------------------------------
	[u m]  = size(U.u);
	[v l]  = size(Y.y);
	n      = length(M.x);

	% get constraints on Cov{e} - C1
	%-------------------------------------------------------------------
	if ~isfield(Y,'Ce')
		Ce = spm_Ce(v*ones(1,l));
	else
		Ce = Y.Ce;
	end

	% Gauss-Newton/Bayesian/EM estimation
	%===================================================================
	[Ep,Cp,Ce] = spm_nlsi_GN(M,U,Y,Ce);

else
	% Use prior expectation to expand around
	%-------------------------------------------------------------------
	Ep         = M.pE;
end

% Bilinear representation
%===========================================================================
[A,B,C,D,L,O] = spm_bireduce(M,Ep);

% stable point ~ x(1e8)
%---------------------------------------------------------------------------
x0            = M.x;
x0            = spm_bilinear(A,B,C,D,x0,1,1e8);


% Volterra kernels
%===========================================================================

% state variables
%---------------------------------------------------------------------------
N             = M.N;
[H0,H1,H2]    = spm_bilinear(A,B,C,D,x0,M.N,M.dt);


% outputs
%---------------------------------------------------------------------------
K0    = feval(M.lx,x0,Ep);
K1    = zeros(N,l,m);
K2    = zeros(N,N,l,m,m);
for i = 1:m
	for p = 1:l
		K1(:,p,i)   = H1(:,:,i)*L(:,p);
	end
end
for i = 1:m
for j = 1:m
	for p = 1:l
		K2(:,:,p,i,j) = H1(:,:,i)*O(:,:,p)*H1(:,:,j)';
		for q = 1:n
			K2(:,:,p,i,j) = K2(:,:,p,i,j) + L(q,p)*H2(:,:,q,i,j);
		end
	end
end
end

% output arguments
%---------------------------------------------------------------------------
if nargin == 3
	varargout = {Ep,Cp,Ce,H0,H1,H2,K0,K1,K2,A,B,C,D};
else
	varargout = {H0,H1,H2,K0,K1,K2,A,B,C,D};
end

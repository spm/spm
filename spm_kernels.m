function [K0,K1,K2,H1] = spm_kernels(M0,M1,L,N,dt)
% returns global Volterra kernels for a MIMO Bilinear system
% FORMAT [K0,K1,K2] = spm_kernels(M0,M1,L,N,dt)
% M0    - (n x n)     df(x(0),0)/dx                    - n states
% M1    - {m}(n x n)  d2f(x(0),0)/dxdu                 - m inputs
% L     - (n x l)     dldx                             - l outputs
% N     - kernel depth       {intervals}
% dt    - interval           {seconds}
%
% Volterra kernels:
%---------------------------------------------------------------------------
% K0    - (1 x l)             = K0(t)         = y(t)
% K1    - (N x l x m)         = K1i(t,s1)     = dy(t)/dui(t - s1)
% K2    - (N x N x l x m x m) = K2ij(t,s1,s2) = d2y(t)/dui(t - s1)duj(t - s2)
%
%___________________________________________________________________________
% Returns Volterra kernels for bilinear systems of the form
%
% dx/dt = f(x,u) = M0*x + M1{1}*x*u1 + ... M1{m}*x*um
%  y(t) = L*x(t)
%
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% Volterra kernels for bilinear systems
%===========================================================================

% parameters
%---------------------------------------------------------------------------
n     = size(M0,1);					% state variables
m     = size(M1,2);					% inputs
l     = size(L ,1);					% ouputs
H1    = zeros(N,n,m);
K1    = zeros(N,l,m);
K2    = zeros(N,N,l,m,m);
M0    = full(M0);


% pre-compute exponentials
%---------------------------------------------------------------------------
e1    = sparse(expm( dt*M0));
e2    = sparse(expm(-dt*M0));
for p = 1:m
	M{1,p} = e1*M1{p}*e2;
end
for i = 2:N
	for p = 1:m
		M{i,p} = e1*M{i - 1,p}*e2;
	end
end

% 0th order kernel
%---------------------------------------------------------------------------
X0    = sparse(1,1,1,n,1);
if nargout > 0
	H0    = e1^N*X0;
	K0    = L*H0;
end


% 1st order kernel
%---------------------------------------------------------------------------
if nargout > 1
	for p = 1:m
	for i = 1:N

		% 1st order kernel
		%-----------------------------------------------------------
		H1(i,:,p) = M{i,p}*H0;
		K1(i,:,p) = H1(i,:,p)*L';
	end
	end
end

% 2nd order kernels
%---------------------------------------------------------------------------
if nargout > 2
	for p = 1:m
	for q = 1:m
	for j = 1:N

		% 2nd order kernel
		%-----------------------------------------------------------
		H  = L*M{j,q}*H1([j:N],:,p)';
		K2(j,[j:N],:,q,p) = H';
		K2([j:N],j,:,p,q) = H';

		% overflow
		%-----------------------------------------------------------
		if max(max(abs(H1([j:N],:,p)))) < 1e-2
			break
		end
	end
	end
	end
end

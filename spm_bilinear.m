function [H0,H1,H2] = spm_bilinear(A,B,C,D,x0,N,dt)
% returns global Volterra kernels for a MIMO Bilinear system
% FORMAT [H0,H1,H2] = spm_bilinear(A,B,C,D,x0,N,dt)
% A     - (n x n)     df(x(0),0)/dx                    - n states
% B     - (n x n x m) d2f(x(0),0)/dxdu                 - m inputs
% C     - (n x m)     df(x(0),0)/du - d2f(x(0),0)/dxdu*x(0)
% D     - (n x 1)     f(x(0).0) - df(x(0),0)/dx*x(0)
% x0    - (n x 1)     x(0)
% N     - kernel depth       {intervals}
% dt    - interval           {seconds}
%
% Volterra kernels:
%
% H0    - (n)       	      = h0(t)       = y(t)
% H1    - (N x n x m)         = h1i(t,s1)    = dy(t)/dui(t - s1)
% H2    - (N x N x n x m x m) = h2ij(t,s1,s2) = d2y(t)/dui(t - s1)duj(t - s2)
%
%___________________________________________________________________________
% Returns Volterra kernels for bilinear systems of the form
%
% dx/dt = f(x,u) = A*x + B1*x*u1 + ... Bm*x*um + C1u1 + ... Cmum + D
%  y(t) = x(t)
%
%---------------------------------------------------------------------------
% %W% Karl Friston %E%


% Volterra kernels for bilinear systems
%===========================================================================

% parameters
%---------------------------------------------------------------------------
n     = size(A,1);					% state variables
m     = size(C,2);					% inputs
t     = N*dt;

% Lie operators
%---------------------------------------------------------------------------
M0    = [0 zeros(1,n); D A];
for i = 1:m
	M1(:,:,i) = [0 zeros(1,n); C(:,i) B(:,:,i)];
end
X0    = [1; x0];

% 0th order kernel
%---------------------------------------------------------------------------
H0    = expm(t*M0)*X0;
if nargout == 1
	H0    = H0([1:n] + 1);
	return
end

% 1st order kernel
%---------------------------------------------------------------------------
H1    = zeros(N,n + 1,m);
for p = 1:m
    for i = 1:N
	u1         = N - i + 1;
	H1(u1,:,p) = expm(u1*dt*M0)*M1(:,:,p)*expm(-u1*M0)*H0;
    end
end
if nargout == 2
	H0    = H0([1:n] + 1);
	H1    = H1(:,[1:n] + 1,:);
	return
end

% 2nd order kernels
%---------------------------------------------------------------------------
H2    = zeros(N,N,n + 1,m,m);
for p = 1:m
for q = 1:m
    for j = 1:N
	u2         = N - j + 1;
        u1         = N - [1:j] + 1;
	H          = expm(u2*dt*M0)*M1(:,:,q)*expm(-u2*dt*M0)*H1(u1,:,p)';
	H2(u2,u1,:,q,p) = H';
	H2(u1,u2,:,p,q) = H';
    end
end
end

% remove kernels associated with the constant state 
%---------------------------------------------------------------------------
H0    = H0([1:n] + 1);
H1    = H1(:,[1:n] + 1,:);
H2    = H2(:,:,[1:n] + 1,:,:);

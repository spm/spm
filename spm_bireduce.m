function [A,B,C,D,L,O] = spm_bireduce(M,P)
% reduction of a fully nonlinear MIMO system to a Bilinear form
% FORMAT [A,B,C,D,L,O] = spm_bireduce(M,P)
%
% M   - model specification structure
% Required fields:
%   M.fx  - dx/dt   = f(x,u,P)  {function string or m-file}
%   M.lx  - y(t)    = l(x,P)    {function string or m-file}
%   M.m   - m inputs
%   M.n   - n states
%   M.l   - l ouputs
%   M.x   - (n x 1) = x(0)
% P     - model parameters
%
% A   - df/dx;
% B   - df/dxu;
% C   - df/du - df/dxu*x0;
% D   - fx0   - df/dx*x0;
% L   - dl/dx
% O   - dl/dxdx
%
%___________________________________________________________________________
% Returns Matrix operators for the Bilinear approximation to the MIMO
% system described by
%
%			dx/dt = f(x,u,P)
% 			 y(t) = l(x,P)
%
% evaluated at x(0) = x and u
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% create inline functions
%---------------------------------------------------------------------------
funx   = fcnchk(M.fx,'x','u','P');
funl   = fcnchk(M.lx,'x','P');

% compute Jacobian and partial derivatives involving x and u
%===========================================================================
m      = M.m;					% m inputs
n      = M.n;					% n states
l      = M.l;					% l ouputs
x      = M.x(:);
dx     = 1e-3;
du     = 1e-3;

% f(x(0),0)
%---------------------------------------------------------------------------
u      = zeros(m,1);
fx0    = feval(funx,x,u,P);

% df(x,0)/du
%---------------------------------------------------------------------------
dfdu   = zeros(n,m);
for  i = 1:m
	ui         = u;
	ui(i)      = ui(i) + du;
	fxu(:,i)   = feval(funx,x,ui,P);
	dfdu(:,i)  = (fxu(:,i) - fx0)/du;
end

% df(x,0)/dx
%---------------------------------------------------------------------------
dfdx   = zeros(n,n);
dfdxu  = zeros(n,n,m);
for  i = 1:n
	xi         = x;
	xi(i)      = xi(i) + dx;
	fxi        = feval(funx,xi,u,P);
	dfdx(:,i)  = (fxi - fx0)/dx;

	% df(x,0)/dxdu
	%-----------------------------------------------------------
	for j = 1:m
		uj           = u;
		uj(j)        = uj(j) + du;
		fxij         = feval(funx,xi,uj,P);
		dfdxj        = (fxij - fxu(:,j))/dx;
		dfdxu(:,i,j) = (dfdxj - dfdx(:,i))/du;
	end
end

% bilinear approximation
%===========================================================================

% matrix operators 
%---------------------------------------------------------------------------
A      = dfdx;
B      = dfdxu;
C      = dfdu;
D      = fx0 - dfdx*x;
for  i = 1:m
	C(:,i) = C(:,i) - dfdxu(:,:,i)*x;
end

% derivatives involving l(x) for output kernels
%===========================================================================
L      = zeros(n,l);
O      = zeros(n,n,l);

if nargout <= 4 | ~l return, end

%    L = dl(x(t))/dx
%---------------------------------------------------------------------------
for  i = 1:n
    xi      = x;
    xi(i)   = xi(i) + dx;
    dldx    = (feval(funl,xi,P) - feval(funl,x,P))/dx;
    L(i,:)  = dldx(:)';

	% O = dl(x(t))/dxdx
	%-------------------------------------------------------------------
	for   j = 1:n
		xj       = x;
		xij      = xi;
		xj(j)    = xj(j)  + dx;
		xij(j)   = xij(j) + dx;
		dldxj    = feval(funl,xij,P);
		dldxj    = (dldxj - feval(funl,xj,P))/dx;
		O(i,j,:) = (dldxj(:) - dldx(:))'/dx;
	end

end

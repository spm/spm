function [y,dy] = spm_int(P,M,U,v)
% integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D;
% FORMAT [y,dy] = spm_int(P,M,U,v)
% P   - model parameters
% M   - model structure
% U   - input structure
% v   - number of sample points [default = 256]
%
% y   - (v x l)  response y = l(x,P)
% dy  - (v x 1)  first temporal derivative dy/dt
%___________________________________________________________________________
% Integrates the bilinear approximation to the MIMO system described by
%  
%    dx/dt = f(x,u,P) = A*x + u*B*x + C*u + D
%    y(t)  = l(x,P)   = L*x; 
%  
% at v sampling points over the input time
%---------------------------------------------------------------------------
% %W% Karl Friston %E%


% embed parameters in matrix operators
%---------------------------------------------------------------------------
if nargin < 4,
	v  = 256;
end

% output nonlinearity
%---------------------------------------------------------------------------
if isfield(M,'lx')
	lx = fcnchk(M.lx,'x','P');
end	

% Bilinear approximation (1st order)
%---------------------------------------------------------------------------
[M0,M1,L]  = spm_bi_reduce(M,P,1);
n          = size(L,2) - 1;			% n states
m          = size(U.u,2);			% m inputs
l          = size(L,1);				% l outputs
u          = size(U.u,1);			% input times

% evaluation time points (when response is sampled or input changes)
%---------------------------------------------------------------------------
s      = [1:v]*u/v;				% output times
t      = find(any(diff(U.u),2))';		% input  times
[T s]  = sort([s t]);				% update (input & ouput) times
dt     = [U.dt*diff(T) 0];			% update intervals
U      = U.u(t + 1,:);

% Integrate
%---------------------------------------------------------------------------
q      = length(dt);
y      = zeros(l,v);
dy     = zeros(l,v);
x      = sparse(1,1,1,n + 1,1);
J      = M0;
if isfield(M,'lx')
	y(:,1) = feval(lx,M.x,P);
else
	y(:,1) = L*x;
end
for  i = 1:q

	% input - update J
	%------------------------------------------------------------------
	if s(i) > v

		u     = U(s(i) - v,:);
		J     = M0;
		for j = 1:m
			J  = J + u(j)*M1{j};
		end

	% output - implement l(x)
	%-------------------------------------------------------------------
	else
		if isfield(M,'lx')
		 	y(:,s(i))  = feval(lx,M.x + x([1:n] + 1),P);
		else
			y(:,s(i))  = L*x;
		end

		if nargout > 1
			dy(:,s(i)) = L*J*x;
		end
	end

	% compute updated states x = expm(J*dt)*x;
        %-------------------------------------------------------------------
        x  = spm_expm(J*dt(i),x);

end
y      = real(y');
dy     = real(dy');

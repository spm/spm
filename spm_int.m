function [y,dy] = spm_int(P,M,U,v,w)
% integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D;
% FORMAT [y,dy] = spm_int(P,M,U,v,w)
% P   - model paramters
% M   - model structure
% U   - input structure
% v   - number of sample points [default = 256]
% w   - delay {secs}            [default = 0  ]
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
	v = 256;
end
if nargin < 5,
	w = 0;
end
n      = M.n;					% n states
m      = M.m;					% m inputs
l      = M.l;					% l outputs
u      = size(U.u,1);				% input times

% output nonlinearity
%---------------------------------------------------------------------------
lx     = fcnchk(M.lx,'x','P');			

% Bilinear approximation (1st order)
%---------------------------------------------------------------------------
[M0,M1,L] = spm_bi_reduce(M,P,1);

% evaluation time points (when response is sampled or input changes)
%---------------------------------------------------------------------------
s      = [1:v]*u/v + w/U.dt;			% output times
t      = find(any(diff(U.u),2))';		% input  times
[T s]  = sort([s t]);				% update (input & ouput) times
dt     = [U.dt*diff(T) 0];			% update intervals
U      = U.u(t + 1,:);

% Integrate
%---------------------------------------------------------------------------
q      = length(dt);
y      = zeros(l,v);
dy     = zeros(l,v);
x      = sparse(1,1,1,length(M0),1);
y(:,1) = feval(lx,M.x,P);
J      = M0;
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
		 y(:,s(i)) = feval(lx,M.x + x([1:n] + 1),P);
		dy(:,s(i)) = L*J*x;
	end

	% compute updated states x = expm(J*dt)*x;
        %-------------------------------------------------------------------
        x  = spm_expm(J*dt(i),x);

end
y      = y';
dy     = dy';

function [y,dy] = spm_int(P,M,U,v)
% integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D;
% FORMAT [y,dy] = spm_int(P,M,U,v)
% P   - model parameters
% M   - model structure
% U   - input structure or matrix
% v   - number of sample points [default = size(U.u,1)]
%
% y   - (v x l)  response y = g(x,u,P)
% dy  - (v x 1)  first temporal derivative dy/dt
%___________________________________________________________________________
% Integrates the bilinear approximation to the MIMO system described by
%  
%    dx/dt = f(x,u,P) = A*x + u*B*x + C*u + D
%    y     = g(x,u,P) = L*x; 
%  
% at v sampling points over the input time
%
% spm_int will also handle static observation models by evaluating
% g(x,u,P)
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% convert U to U.u if necessary
%---------------------------------------------------------------------------
if ~isstruct(U),
	U.u  = U;
end
try
	U.dt;
catch
	U.dt = 0;
end


% add [0] states if not specified
%---------------------------------------------------------------------------
if ~isfield(M,'f')
	M.f = inline('sparse(0,1)','x','u','P');
	M.n = 0;
	M.x = sparse(0,0);
end

% number of times to sample
%---------------------------------------------------------------------------
if nargin < 4,
	v  = size(U.u,1);
else
	v  = min([v size(U.u,1)]);
end

% output nonlinearity, if specified
%---------------------------------------------------------------------------
if isfield(M,'g')
	g  = fcnchk(M.g,'x','u','P');
end	

% Bilinear approximation (1st order)
%---------------------------------------------------------------------------
[M0,M1,L]  = spm_bi_reduce(M,P);
n          = size(L,2) - 1;			% n states
m          = size(U.u,2);			% m inputs
l          = size(L,1);				% l outputs
u          = size(U.u,1);			% input times

% evaluation time points (when response is sampled or input changes)
%---------------------------------------------------------------------------
s      = ceil([1:v]*u/v);			% output times
t      = [1 (1 + find(any(diff(U.u),2))')];	% input  times
[T s]  = sort([s t]);				% update (input & ouput) times
dt     = [U.dt*diff(T) 0];			% update intervals

% Integrate
%---------------------------------------------------------------------------
y      = zeros(l,v);
dy     = zeros(l,v);
x      = sparse(1,1,1,n + 1,1);
J      = M0;
for  i = 1:length(T)

	% input
	%------------------------------------------------------------------
	u     = U.u(T(i),:);

	% change in input - update J
	%------------------------------------------------------------------
	if s(i) > v

		J     = M0;
		for j = 1:m
			J  = J + u(j)*M1{j};
		end

	% output sampled - implement l(x)
	%-------------------------------------------------------------------
	else
		if isfield(M,'g')
		 	y(:,s(i))  = feval(g,M.x + x([1:n] + 1),u,P);
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

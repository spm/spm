function [y] = spm_nlsi_int(P,M,U,v,w)
% integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D;
% FORMAT [y] = spm_nlsi_int(P,M,U,v,w)
% P   - model  paramters
% M   - model  structure
% U   - input  structure
% v   - number of sample points [default = 256]
% w   - delay {secs}            [default = 0  ]
%
% y     - (v x l)  response y = l(x,P)
%___________________________________________________________________________
% Integrates the bilinear approximation to the MIMO system described by
%  
%    dx/dt = f(x,u,P) = A*x + u*B*x + C*u + D
%    y(t)  = l(x,P)   = x*L + x'*O*x/2 
%  
% at v sampling points over the input time
%---------------------------------------------------------------------------
% %W% Karl Friston %E%


% embed parameters in matrix operators
%---------------------------------------------------------------------------
if nargin < 4,
	v    = 256;
end
if nargin < 5,
	w    = 0;
end
n            = M.n;					% n states
m            = M.m;					% m inputs
l            = M.l;					% l outputs
u            = size(U.u,1);				% input times

% Bilinear approximation
%---------------------------------------------------------------------------
[A,B,C,D]    = spm_bireduce(M,P);

% stable point x(0)
%---------------------------------------------------------------------------
x0           = spm_bilinear(A,B,C,D,M.x,1,v*U.dt);
y0           = feval(M.lx,x0,P);

% make Lie matrices sparse
%---------------------------------------------------------------------------
A            = sparse(A);
C            = sparse(C);
D            = sparse(D);
for i = 1:m
	E{i} = sparse(B(:,:,i));
end

% evaluation time points (when response is sampled or input changes)
%---------------------------------------------------------------------------
s      = [1:v]*u/v + w/U.dt;			% output times
t      = find(any(diff(U.u),2))';		% input  times
[T s]  = sort([s t]);				% update times
dt     = U.dt*diff(T);				% update intervals
U      = U.u(t + 1,:);		

% Integrate
%---------------------------------------------------------------------------
q      = length(dt);
y      = zeros(v,l);
x      = x0;
y(1,:) = y0;
J      = A;
K      = D;
for  i = 1:q

	% input - update J and K
	%------------------------------------------------------------------
	if s(i) > v

		u     = U(s(i) - v,:);
		J     = A;
		for j = 1:m
			J  = J + u(j)*E{j};
		end
		K     = C*u' + D;

	% output - implement l(x)
	%-------------------------------------------------------------------
	else
		y(s(i),:) = feval(M.lx,x,P);
	end

        % compute dx = (expm(J*dt) - eye(n))*inv(J)*fx;
        %-------------------------------------------------------------------
        fx    = (J*x + K)*dt(i);
        dx    = fx;
        for j = 1:32

                fx = J*fx*dt(i)/(j + 1);
                dx = dx + fx;

		% break if converging
		%-----------------------------------------------------------
		if fx'*fx < 1e-16
			break
		end
        end

	% update
	%-------------------------------------------------------------------
        x     = x + dx;

end

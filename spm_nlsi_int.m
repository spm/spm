function [y] = spm_nlsi_int(P,M,U,v)
% integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D;
% FORMAT [y] = spm_nlsi_int(P,M,U,v)
% P   - model  paramters
% M   - model  structure
% U   - input  structure
% v   - number of sample points [default = 256]
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
n            = M.n;					% n states
m            = M.m;					% m inputs
l            = M.l;					% l outputs
u            = size(U.u,1);				% input times

% Bilinear approximation
%---------------------------------------------------------------------------
[A,B,C,D]    = spm_bireduce(M,P);

% stable point x(0)
%---------------------------------------------------------------------------
x0           = spm_bilinear(A,B,C,D,M.x,1,1e8);		% x(1e8)
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
s      = ceil([1:v]*u/v);
s      = sparse(s,1,1,u,1);			% observation times
t      = ~[~s & [0; all(~diff(U.u),2)]];	% time points
dt     = U.dt*diff(find(t));			% time intervals
s      = s(t);					% sample times
u      = U.u(t,:);

% Preliminary computations
%---------------------------------------------------------------------------
C      = C*u';
for  i = 1:n
	C(i,:) = C(i,:) + D(i);
end

% Integrate
%---------------------------------------------------------------------------
q      = length(dt);
y      = zeros(q + 1,l);
x      = x0;
y(1,:) = y0;
for  i = 1:q

	% f(x)
        %-------------------------------------------------------------------
        J     = A;
        for j = 1:m
                J  = J + u(i,j)*E{j};
        end

        % compute dx = (expm(J*dt) - eye(n))*inv(J)*fx;
        %-------------------------------------------------------------------
        fx    = (J*x + C(:,i))*dt(i);
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
        x     = x + dx;

	% implement l(x)
	%-------------------------------------------------------------------
	if s(i + 1)
		y(i + 1,:) = feval(M.lx,x,P);
	end 

end

% Subsample 
%---------------------------------------------------------------------------
y     = y(find(s),:);

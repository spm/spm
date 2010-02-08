
% Notes for Laplace scheme and dynamic updating of parameters
%==========================================================================
n  = 128;                       % number of samples (time bins)
P  = 4;                         % true parameter
k  = n;                         % precision on fluctuations
pp = 2;                         % prior on parameter
m  = 8;                         % number of data
iS = eye(m,m)*2;                % error precision
s  = sqrtm(inv(iS));
 
p  = [0;0];                     % initial parameter estimates
for i = 1:n
    
    x      = 4 + randn(m,1);          % exogenous input
    z      = s*randn(m,1);            % error
    y      = P*x + z;                 % response
    e      = y - p(1)*x;              % prediction error
    Lp     = -e'*iS*x + pp*p(1);
    Lpp    = x'*iS*x  + pp;
    f      = [p(2); (-Lp -k*p(2))];
    dfdx   = [0 1;   -Lpp -k];
    p      = p + spm_dx(dfdx,f,1);
    X(:,i) = p;
    LP(i)  = Lp;
end
 
% results
%--------------------------------------------------------------------------
subplot(2,1,1)
plot(1:n,X)
title('generalised parameter estimates','FontSize',16)
 
subplot(2,1,2)
plot(1:n,LP)
title('energy gradient','FontSize',16)

 
return

% Notes
%==========================================================================
M.f = inline('[x(2); (u - K(1)*x(2))]','x','u','K','M');
M.m = 1;
M.n = 2;
M.l = 2;
M.x = [0;0];
M.u = 0;

N   = 128;
dt  = 1/64;
K   = 4;
[K0,K1,K2] = spm_kernels(M,K,N,dt);

subplot(2,1,1)
plot([1:N]*dt,K1);
xlabel('time (s)')
title('Kernels','FontSize',16)
legend('drive','trace')
axis square


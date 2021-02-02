function FEP_lorenz_surprise
% This demo computes the cost-function (negative reward) for a Lorenz
% system; to show cost can be computed easily from value (negative 
% surprise or sojourn time). However, value is not a Lyapunov function
% because the flow is not curl-free (i.e., is not irrotational).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: FEP_lorenz_surprise.m 8046 2021-02-02 18:48:05Z karl $


%% generative model - linearised Lorentz system for doublechecking
%==========================================================================                       % switch for demo
% f   = @(x,v,P,M) [-P(1) 0 0; 0 -P(2) 0; 0 0 -P(3)]*x/64;
% P   = [1; 1; 8];
% x0  = [0; 0; 0];
% W   = diag([1; 1; 1/4]);
 
% dynamics and parameters of a Lorentz system (with Jacobian)
%==========================================================================
% dxdt = f(x) + w
%--------------------------------------------------------------------------
f    = @(x,v,P,M) [P(1)*x(2) - P(1)*x(1);
                   P(3)*x(1) - x(2) - x(1)*x(3);
                   P(2)*x(3) + x(1)*x(2)]/64;
J    = @(x,v,P,M) [[     -P(1),P(1),     0];
                   [P(3) - x(3), -1, -x(1)];
                   [     x(2), x(1),  P(2)]]/64;
P    = [10; -8/3; 32];                 % parameters
x0   = [1; 1; 32];                     % initial states
W    = diag([1 1 1]/64);               % precision of random fluctuations


% coupled oscillators
%==========================================================================
% dxdt = f(x) + w
%--------------------------------------------------------------------------
% f    = @(x,v,P,M) [
%     P(1)*x(2) - P(1)*(x(1)*(1 - 1/8) + x(4)/8);
%     P(3)*x(1) - x(2) - x(1)*x(3);
%     P(2)*x(3) + x(1)*x(2);
%     P(1)*x(5) - P(1)*(x(4)*(1 - 1/8) + x(1)/8);
%     P(3)*x(4) - x(5) - x(4)*x(6);
%     P(2)*x(6) + x(4)*x(5)]/64;
% 
% J    = @(x,v,P,M) [
%     [P(1)*(c - 1), P(1),   0,  -P(1)/8,  0,   0]
%     [  P(3) - x(3), -1, -x(1),       0,  0,   0]
%     [  x(2), x(1),  P(2),            0,  0,   0]
%     [ -P(1)/8,  0,   0, P(1)*(c - 1), P(1),   0]
%     [       0,  0,   0,    P(3) - x(6), -1, -x(4)]
%     [       0,  0,   0,         x(5), x(4),  P(2)]]/64;
% 
% P    = [10; -8/3; 32];                 % parameters
% x0   = [1; 1; 32; 0; 0; 16];           % initial states
% W    = diag([1 1 1 1 1 1]/16);         % precision of random fluctuations


% state space model
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;
M.W  = W;
M.K  = 3;

% solution of SDE
% -------------------------------------------------------------------------
spm_figure('GetWin','NESS (1)'); clf
T    = 2^10;
U.u  = sparse(T,M(1).m);
U.dt = 1;
t    = spm_int_ode(M.pE,M,U);
plot(t)

% state space
%--------------------------------------------------------------------------
N    = 16;                                     % number of bins
x{1} = linspace(-32,32,N);
x{2} = linspace(-32,32,N);
x{3} = linspace( 16,48,N);

% x{4} = linspace(-32,32,N);
% x{5} = linspace(-32,32,N);
% x{6} = linspace( 16,48,N);

% Fokker-Planck operator and equilibrium density
%==========================================================================
M.X             = spm_ndgrid(x);               % sample points
[q0,X,F,f,NESS] = spm_ness_GN(M);              % NESS density
nx              = size(q0);                    % coarse graining            
n               = numel(nx);                   % dimensionality

% solve for a trajectory 
%--------------------------------------------------------------------------
T     = 1024;
s     = x0;
dt    = 1/2;
t     = zeros(n,T);
for i = 1:T
    d      = bsxfun(@minus,X,s');
    d      = sum(d.^2,2);
    p      = spm_softmax(-d/16);
    ds     = F'*p;
    s      = s + ds*dt; % + sqrtm(inv(W))*randn(n,1)*dt;
    t(:,i) = s;
end

% solve for a trajectory 
%--------------------------------------------------------------------------
Ep    = NESS.Ep;
s     = x0;
for i = 1:T
    M.X    = s';
    ds     = spm_NESS_gen(Ep,M);
    s      = s + ds'*dt; % + sqrtm(inv(W))*randn(n,1)*dt;
    t(:,i) = s;
end



% illustrate convergence to nonequilibrium steady-state
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (1)'); clf
subplot(3,2,3)
plot(t')
axis square
title('trajectory','Fontsize',16)
xlabel('time'), ylabel('states')

subplot(3,2,4)
plot3(t(1,:),t(2,:),t(3,:))
% hold on, plot3(t(4,:),t(5,:),t(6,:)), hold off
axis square, title('trajectory','Fontsize',16)


% marginal nonequilibrium steady-state
%--------------------------------------------------------------------------
% p0 = spm_marginal(q0);


% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,1)
imagesc(x{3},x{2},squeeze(sum(q0,1)))
hold on, plot(t(3,:),t(2,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,2)
imagesc(x{3},x{1},squeeze(sum(q0,2)))
hold on, plot(t(3,:),t(1,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,3)
imagesc(x{2},x{1},squeeze(sum(q0,3)))
hold on, plot(t(2,:),t(1,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

subplot(3,1,3)
for i = 1:3
    plot(f(i,:)',F(:,i),'.','MarkerSize',1), hold on
end, hold off
axis square xy
title('flows','Fontsize',16)
xlabel('approximate flow'), ylabel('true flow')


%% illustrate conditional dependencies
%==========================================================================
spm_figure('GetWin','NESS (2)'); clf

% covariance
%--------------------------------------------------------------------------
m    = q0(:)'*X;
C    = full(X'*bsxfun(@times,q0(:),X) - m'*m);

% log normal of Jacobian and hessian
%--------------------------------------------------------------------------
subplot(3,3,1)
imagesc(log(NESS.J2 + exp(-16))), axis square
title('log |Jacobian|','Fontsize',12)

subplot(3,3,2)
imagesc(log(inv(C).^2 + exp(-16))), axis square
title('log |Precision|','Fontsize',12)

subplot(3,3,3)
imagesc(log(NESS.H2 + exp(-16))), axis square
title('log |Hessian|','Fontsize',12)

% illustrate entities - conditioning on blanket state (1st state)
%--------------------------------------------------------------------------
N     = 64;
x{1} = linspace(-32,32,N);
x{2} = linspace(-32,32,N);
x{3} = linspace( 16,48,N);

nx    = [N N N];
p0    = spm_softmax(spm_polymtx(x,M.K)*NESS.bS);
p0    = reshape(p0,[N N N]);

for i = 1:8
    subplot(6,4,8 + i)
    j = ceil(i*(nx(1)/8));
    imagesc(x{3},x{2},squeeze(p0(j,:,:)).^8)
    axis square xy
    title(sprintf('Conditional x = %.2g',x{1}(j)),'Fontsize',12)
end

return


%% Notes for functional form of polynomial flow
%--------------------------------------------------------------------------
M.X         = M.x(:)';
[F,Q,S,L,H] = spm_NESS_gen(NESS.Ep,M)

clear Ep
syms x1 x2 x3 W 'real'
M.X         = [x1,x2,x3];
M.W         = diag([W,W,W]);

for i = 1:numel(NESS.Ep.Sp)
    if abs(NESS.Ep.Sp(i)) > exp(-8)
        Ep.Sp(i,1) = sym(sprintf('S%i',i),'real');
    else
        Ep.Sp(i,1) = 0;
    end
end
for i = 1:numel(NESS.Ep.Qp)
    if abs(NESS.Ep.Qp(i)) > exp(-1)
        Ep.Qp(i,1) = sym(sprintf('Q%i',i),'real');
    else
        Ep.Qp(i,1) = 0;
    end
end

U           = spm_ness_U(M);
[F,Q,S,L,H] = spm_NESS_gen(Ep,M)
Q           = reshape(cat(1,Q{:}),n,n)
f           = real(U.f)



%% illustrate extrinsic information geometry
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (3)'); clf

q     = zeros(N,T);
for i = 1:T
    
    % conditional (marginal) probability over second state, given the first
    %----------------------------------------------------------------------
    [m,j]  = min(abs(x{1} - t(1,i)));
    q2     = squeeze(spm_sum(p0(j,:,:),3));
    q2     = q2(:)/sum(q2);
    q(:,i) = q2;
    
    q3     = squeeze(spm_sum(p0(j,:,:),2));
    q3     = q3(:)/sum(q3);
    
    % conditional manifold (conditional expectations of second and third)
    %----------------------------------------------------------------------
    m2(i)   = x{2}*q2;
    m3(i)   = x{3}*q3;
    
end

subplot(4,2,1)
plot(t'), spm_axis tight

subplot(4,2,2)
plot(m2,m3,'.')

subplot(4,1,2)
imagesc(1:T,x{2},(1 - q)), hold on
plot(t(2,:),'r'), plot(m2,'r:'), axis xy
title('Conditional','Fontsize',14)

subplot(4,1,3)
plot(m3,'r:')
title('Conditional','Fontsize',14)

subplot(4,1,4)
plot(t(3,:),'r')
title('Conditional','Fontsize',14)


return



%% Notes (symbolic maths or Lorenz system)
%==========================================================================
syms x1 x2 x3 x4 x5 x6 P1 P2 P3 f J x P c

x    = [x1;x2;x3];
v    = [x4;x5;x6];
P    = [P1;P2;P3];

f    = symfun([P(1)*x(2) - P(1)*x(1);
               P(3)*x(1) - x(2) - x(1)*x(3);
               P(2)*x(3) + x(1)*x(2)]/64,x);
    
J    = [diff(f,x1) diff(f,x2) diff(f,x3)]
K    = [diff(J,x1) diff(J,x2) diff(J,x3)]
f    = f(0,0,0) + J*x + (1/2)*K*kron(x,x)

% generalised coordinates
%--------------------------------------------------------------------------
ff   = symfun([v; J*v + (1/2)*K*kron(x,v) + (1/2)*K*kron(v,x)],[x;v])
JJ   = [diff(ff,x1) diff(ff,x2) diff(ff,x3) diff(ff,x4) diff(ff,x5) diff(ff,x6)]

% coupled oscillators
%--------------------------------------------------------------------------
x    = [x1;x2;x3;x4;x5;x6];
ff   = symfun([P(1)*x(2) - P(1)*(x(1)*(1 - c) + x(4)*c);
               P(3)*x(1) - x(2) - x(1)*x(3);
               P(2)*x(3) + x(1)*x(2);
               P(1)*x(5) - P(1)*(x(4)*(1 - c) + x(1)*c);
               P(3)*x(4) - x(5) - x(4)*x(6);
               P(2)*x(6) + x(4)*x(5)]/64,x);


J   = [diff(ff,x1) diff(ff,x2) diff(ff,x3) diff(ff,x4) diff(ff,x5) diff(ff,x6)]


%% generalised coordinates of motion
%==========================================================================
% dxdt = x' + w
% dx'dt = J(x)*x' + ... + w'
% J     = dx'dx
% |dxdt | = |x'     | + |w |
% |dx'dt| = |J(x)*x' + ... | + |w'|
% J'      = |0 ...|
%--------------------------------------------------------------------------
f    = @(x,v,P,M) [                               x(4)
                                                  x(5)
                                                  x(6)
                             (P(1)*x(5))/64 - (P(1)*x(4))/64
x(4)*(P(3)/64 - x(3)/64) - (x(1)*x(6))/32 - (x(3)*x(4))/64 - x(5)/64
                (P(2)*x(6))/64 + (x(1)*x(5))/32 + (x(2)*x(4))/32];

J    = @(x,v,P,M) [
[     0,     0,      0,             1,     0,      0]
[     0,     0,      0,             0,     1,      0]
[     0,     0,      0,             0,     0,      1]
[     0,     0,      0,        -P1/64, P1/64,      0]
[-x6/32,     0, -x4/32, P3/64 - x3/32, -1/64, -x1/32]
[ x5/32, x4/32,      0,         x2/32, x1/32,  P2/64]];

P    = [10; -8/3; 32];                 % parameters
x0   = [1; 1; 32; 0; 0; -1];           % initial states
W    = diag([1 1 1 1 1 1]/64);         % precision of random fluctuations


function FEP_lorenz_surprise
% This demo provides an elementary characterisation of stochastic chaos
% using the Lorenz system. Effectively, it uses iterated least-squares to
% solve the Helmholtz decomposition of nonequilibrium steady-state flow
% (i.e., the solution to the Fokker Planck equation) using the Lorentz
% system as an example. This furnishes a generative model for stochastic
% chaos in terms of the underlying potential (log nonequilibrium
% steady-state density) and flow operator, with symmetric and antisymmetric
% (skew symmetric) components. The latter (solenoidal) part of the flow
% operator breaks detailed balance and renders the steady-state a
% nonequilibrium steady-state.
%
% In virtue of using a polynomial expansion for the nonequilibrium
% potential (i.e., surprisal or self information) one can approximate the
% expected flow with a second order polynomial. This can be regarded as a
% Laplace approximation to the nonequilibrium steady-state density. In this
% example, the bivariate distribution of any two states, conditioned upon
% the remainder is Gaussian. Further constraints can be used to specify the
% stochastic chaos as (state dependent) solenoidal flow around a
% multivariate Gaussian, which might be a reasonable approximation in the
% setting of high amplitude random fluctuations.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: FEP_lorenz_surprise.m 8074 2021-03-02 20:58:49Z karl $


%% generative model - linearised Lorenz system for doublechecking
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
x0   = [1; 1; 24];                     % initial states
W    = diag([1 1 1]/16);               % precision of random fluctuations

% state-space model (for SPM integrators)
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;
M.W  = W;

% solution of underlying ordinary differential equation
% -------------------------------------------------------------------------
spm_figure('GetWin','NESS (1)'); clf
T    = 2^10;
U.u  = sparse(T,M(1).m);
U.dt = 1;
t    = spm_int_ode(M.pE,M,U);
plot(t)

% state-space four (Laplace) solution 
%--------------------------------------------------------------------------
N    = 5;                                     % number of bins
d    = 16;
x{1} = linspace(-d,d,N);
x{2} = linspace(-d,d,N);
x{3} = linspace(28 - d,28 + d,N);

% Fokker-Planck operator and equilibrium density
%==========================================================================
[q0,X,F,f,NESS] = spm_ness_hd(M,x);            % NESS density
nx              = size(q0);                    % coarse graining            
n               = numel(nx);                   % dimensionality

% solve for a trajectory using a polynomial approximation to flow 
%--------------------------------------------------------------------------
T     = 1024;
t     = zeros(n,T);
dt    = 1/4;
s     = x0;
Ep    = NESS.Ep;
for i = 1:T
    M.X    = s';
    ds     = spm_NESS_gen(Ep,M);
    s      = s + ds'*dt; 
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
axis square, title('trajectory','Fontsize',16)

% marginal nonequilibrium steady-state
%--------------------------------------------------------------------------
N    = 32;
x{1} = linspace(-32,32,N);
x{2} = linspace(-32,32,N);
x{3} = linspace(  0,64,N);
p0   = spm_softmax(spm_polymtx(x)*NESS.Ep.Sp);
p0   = reshape(p0,[N N N]);

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,1)
imagesc(x{3},x{2},squeeze(sum(p0,1)))
hold on, plot(t(3,:),t(2,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,2)
imagesc(x{3},x{1},squeeze(sum(p0,2)))
hold on, plot(t(3,:),t(1,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,3)
imagesc(x{2},x{1},squeeze(sum(p0,3)))
hold on, plot(t(2,:),t(1,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

subplot(3,1,3)
for i = 1:n
    plot(f(i,:)',F(:,i),'.','MarkerSize',1), hold on
end, hold off
axis square xy
title('Expected flows','Fontsize',16)
xlabel('approximate flow'), ylabel('true flow')


%% illustrate conditional dependencies
%==========================================================================
spm_figure('GetWin','NESS (2)'); clf

% covariance based upon a Laplace approximation to surprisal
%--------------------------------------------------------------------------
% [m,C] = spm_ness_cond(n,3,NESS.Ep.Sp);
m     = q0(:)'*X;
C     = full(X'*bsxfun(@times,q0(:),X) - m'*m);

% log norms of (generally state-dependent) Jacobian and Hessian
%--------------------------------------------------------------------------
subplot(3,3,1)
imagesc(-log(NESS.J2 + exp(-16))), axis square
title('log |Jacobian|','Fontsize',12)

subplot(3,3,2)
imagesc(-log(C.^2 + exp(-16))), axis square
title('log |Covariance|','Fontsize',12)

subplot(3,3,3)
imagesc(-log(NESS.H2 + exp(-16))), axis square
title('log |Hessian|','Fontsize',12)

% illustrate conditional densities
%==========================================================================
N    = 64;
x{1} = linspace(-32,32,N);
x{2} = linspace(-32,32,N);
x{3} = linspace( 16,48,N);

nx  = [N N N];
p0  = spm_softmax(spm_polymtx(x)*NESS.Ep.Sp);
p0  = reshape(p0,nx);


% illustrate marginal densities
%--------------------------------------------------------------------------
s     = round(linspace(-24,24,5));
for i = 1:numel(s)
    subplot(6,numel(s),numel(s)*2 + i)
    [m,j] = min(abs(x{1} - s(i)));
    imagesc(x{3},x{2},squeeze(p0(j,:,:)).^2)
    axis square xy, axis off
    title('2nd and third','Fontsize',12)
end
s     = round(linspace(-24,24,5));
for i = 1:numel(s)
    subplot(6,numel(s),numel(s)*3 + i)
    [m,j] = min(abs(x{2} - s(i)));
    imagesc(x{3},x{1},squeeze(p0(:,j,:)).^2)
    axis square xy, axis off
    title('first and third','Fontsize',12)
end
s     = round(linspace(16,48,5));
for i = 1:numel(s)
    subplot(6,numel(s),numel(s)*4 + i)
    [m,j] = min(abs(x{2} - s(i)));
    imagesc(x{2},x{1},squeeze(p0(:,:,j)).^2)
    axis square xy, axis off
    title('first and 2nd','Fontsize',12)
end

return

%% functional form of polynomial flow
%--------------------------------------------------------------------------
clear Ep
syms  W 'real'
syms  x [1 n] 'real'

% replace sample points with symbolic variables
%--------------------------------------------------------------------------
M.X   = x;
M.W   = diag(kron(ones(n,1),W));

% set small polynomial coefficients to zero
%--------------------------------------------------------------------------
ep    = NESS.Ep;
for i = 1:numel(NESS.Ep.Sp)
    Ep.Sp(i,1) = sym(sprintf('S%i',i),'real');
    if abs(NESS.Ep.Sp(i)) < 1e-4
        Ep.Sp(i,1) = 0;
        ep.Sp(i,1) = 0;
    end
end
for i = 1:numel(NESS.Ep.Qp)
    Ep.Qp(i,1) = sym(sprintf('Q%i',i),'real');
    if abs(NESS.Ep.Qp(i)) < 1e-1
        Ep.Qp(i,1) = 0;
        ep.Qp(i,1) = 0;
    end
end

% evaluate flow, flow operators and Hessians for display
%--------------------------------------------------------------------------
sympref('FloatingPointOutput',0);
sympref('MatrixWithSquareBrackets',true);

U           = spm_ness_U(M);
[F,Q,S,L,H] = spm_NESS_gen(Ep,M);
F = F'
Q = reshape(cat(1,Q{:}),n,n)
S = S
L = L'
H = H
f = real(U.f)

% latex(H)
% sympref('FloatingPointOutput',true);
% sympref('MatrixWithSquareBrackets',true);
% 
% [F,Q,S,L,H] = spm_NESS_gen(ep,M);
% F = F'
% Q = reshape(cat(1,Q{:}),n,n)
% S = S
% L = L'
% H = H
% f = real(U.f)


%% illustrate extrinsic information geometry
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (3)'); clf

% get marginals over external and internal states
%--------------------------------------------------------------------------
q     = zeros(N,T);
b     = 2;                                    % blanket state
h     = 1;                                    % external (hidden) states
m     = 3;                                    % internal states
for i = 1:T
    
    % get marginals over external and internal states
    %----------------------------------------------------------------------
    s       = x;
    s{b}    = t(b,i);
    s{m}    = 0;
    pH(:,i) = spm_softmax(spm_polymtx(s)*NESS.Ep.Sp);
    
    s       = x;
    s{b}    = t(b,i);
    s{h}    = 0;
    pI(:,i) = spm_softmax(spm_polymtx(s)*NESS.Ep.Sp);
    
    % manifold (conditional expectations)
    %-------------------------------------------====-----------------------
    eH(i)   = x{h}*pH(:,i);
    eI(i)   = x{m}*pI(:,i);
    
end

subplot(4,2,2)
plot(eI,eH,'.')

subplot(4,2,1)
plot(eI,t(m,:),'.')

subplot(4,1,2), hold off
imagesc(1:T,x{h},(1 - pH).^2), hold on
plot(t(h,:),'r'), plot(eH,'w'), axis xy
title('Conditional','Fontsize',14)

subplot(4,1,3)
plot(eI,'k')
title('Conditional','Fontsize',14)

return



%% coupled oscillators
%--------------------------------------------------------------------------
% The following notes reproduce the above analysis but for two coupled
% Lorenz systems. Plausible solutions require extra constraints on the
% polynomial approximation to the nonequilibrium steady-state potential.
% These can be implemented by editing the appropriate subroutines to
% enforce a multivariate Gaussian approximation – and eliminate
% second-order terms that are precluded by the Jacobian (i.e., dynamical
% coupling)
%==========================================================================
% dxdt = f(x) + w
%--------------------------------------------------------------------------
f    = @(x,v,P,M) [
    P(1)*x(2) - P(1)*(x(1)*(1 - P(4)) + x(4)*P(4));
    P(3)*x(1) - x(2) - x(1)*x(3);
    P(2)*x(3) + x(1)*x(2);
    P(1)*x(5) - P(1)*(x(4)*(1 - P(4)) + x(1)*P(4));
    P(3)*x(4) - x(5) - x(4)*x(6);
    P(2)*x(6) + x(4)*x(5)]/64;

J    = @(x,v,P,M) [
    [P(1)*(P(4) - 1), P(1),   0,  -P(1)*P(4),  0,   0]
    [  P(3) - x(3), -1, -x(1),             0,  0,   0]
    [  x(2), x(1),  P(2),                  0,  0,   0]
    [ -P(1)*P(4),    0,   0, P(1)*(P(4) - 1), P(1), 0]
    [       0,  0,   0,   P(3) - x(6),      -1, -x(4)]
    [       0,  0,   0,          x(5), x(4),     P(2)]]/64;

P    = [10; -8/3; 32; 1/16];            % parameters
x0   = [1; 1; 32; 1; 1; 28];            % initial states
W    = diag([1 1 1 1 1 1]/8);           % precision of random fluctuations

% state space model
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;
M.W  = W;

% solution of SDE
% -------------------------------------------------------------------------
spm_figure('GetWin','NESS (3)'); clf
T    = 2^10;
U.u  = sparse(T,M(1).m);
U.dt = 1;
t    = spm_int_ode(M.pE,M,U);
plot(t)

% state space
%--------------------------------------------------------------------------
N    = 5;                                     % number of bins
d    = 16;
x{1} = linspace(-d,d,N);
x{2} = linspace(-d,d,N);
x{3} = linspace(28 - d,28 + d,N);
x{4} = linspace(-d,d,N);
x{5} = linspace(-d,d,N);
x{6} = linspace(28 - d,28 + d,N);

% Fokker-Planck operator and equilibrium density
%==========================================================================
[q0,X,F,f,NESS] = spm_ness_hd(M,x);            % NESS density
nx              = size(q0);                    % coarse graining            
n               = numel(nx);                   % dimensionality

% solve for a trajectory 
%--------------------------------------------------------------------------
T     = 1024;
t     = zeros(n,T);
dt    = 1/2;
s     = x0;
Ep    = NESS.Ep;
for i = 1:T
    M.X    = s';
    ds     = spm_NESS_gen(Ep,M);
    s      = s + ds'*dt; % + sqrtm(inv(W))*randn(n,1)*dt;
    t(:,i) = s;
end
clear s

% illustrate convergence to nonequilibrium steady-state
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (3)'); clf
subplot(3,2,3)
plot(t')
axis square
title('trajectory','Fontsize',16)
xlabel('time'), ylabel('states')

subplot(3,2,4)
plot3(t(1,:),t(2,:),t(3,:))
hold on, plot3(t(4,:),t(5,:),t(6,:)), hold off
axis square, title('trajectory','Fontsize',16)


% marginal nonequilibrium steady-state
%--------------------------------------------------------------------------
N    = 32;
s{1} = linspace(-32,32,N);
s{2} = linspace(-32,32,N);
s{3} = linspace( 16,64,N);
s{4} = 0;
s{5} = 0;
s{6} = 0;
p0   = spm_softmax(spm_polymtx(s,M.K)*NESS.Ep.Sp);
p0   = reshape(p0,[N N N]);

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,1)
imagesc(x{3},x{2},squeeze(sum(p0,1)))
hold on, plot(t(3,:),t(2,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,2)
imagesc(x{3},x{1},squeeze(sum(p0,2)))
hold on, plot(t(3,:),t(1,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,3)
imagesc(x{2},x{1},squeeze(sum(p0,3)))
hold on, plot(t(2,:),t(1,:),'r.'), hold off
axis square xy
title('NESS density','Fontsize',16)

subplot(3,1,3)
for i = 1:n
    plot(f(i,:)',F(:,i),'.','MarkerSize',1), hold on
end, hold off
axis square xy
title('flows','Fontsize',16)
xlabel('approximate flow'), ylabel('true flow')


%% illustrate conditional dependencies
%==========================================================================
spm_figure('GetWin','NESS (4)'); clf

% covariance
%--------------------------------------------------------------------------
[m,C] = spm_ness_cond(n,M.K,NESS.Ep.Sp);

% log normal of Jacobian and hessian
%--------------------------------------------------------------------------
subplot(3,3,1)
imagesc(-log(NESS.J2 + exp(-16))), axis square
title('log |Jacobian|','Fontsize',12)

subplot(3,3,2)
imagesc(-log(C.^2 + exp(-16))), axis square
title('log |Covariance|','Fontsize',12)

subplot(3,3,3)
imagesc(-log(NESS.H2 + exp(-16))), axis square
title('log |Hessian|','Fontsize',12)


%% Notes for functional form of polynomial flow
%--------------------------------------------------------------------------
% M.X         = M.x(:)';
% [F,Q,S,L,H] = spm_NESS_gen(NESS.Ep,M)

clear Ep
syms x1 x2 x3 x4 x5 x6 W 'real'

M.X   = [x1,x2,x3,x4,x5,x6];
M.W   = diag([W,W,W,W,W,W]);

for i = 1:numel(NESS.Ep.Sp)
    Ep.Sp(i,1) = sym(sprintf('S%i',i),'real');
    if abs(NESS.Ep.Sp(i)) < exp(-16)
        Ep.Sp(i,1) = 0;
    end
end
for i = 1:numel(NESS.Ep.Qp)
    Ep.Qp(i,1) = sym(sprintf('Q%i',i),'real');
    if abs(NESS.Ep.Qp(i)) < exp(-4)
        Ep.Qp(i,1) = 0;
    end
end

U           = spm_ness_U(M);
[F,Q,S,L,H] = spm_NESS_gen(Ep,M)
Q           = reshape(cat(1,Q{:}),n,n)
f           = real(U.f)



%% illustrate extrinsic information geometry
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (5)'); clf

b     = 4;                                    % blanket state
h     = 1:3;                                  % external (hidden) states
m     = 5:6;                                  % internal states

% get conditional marginals over external and internal states
%--------------------------------------------------------------------------
for i = 1:T
    [m,C]  = spm_ness_cond(n,M.K,NESS.Ep.Sp,b,t(b,i));
    E(:,i) = m;
    V(:,i) = diag(C);
end

subplot(4,2,1)
plot(t'), spm_axis tight


spm_plot_ci(E(1,:),V(1,:)')

return

%% illustrate extrinsic information geometry
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (5)'); clf

q     = zeros(N,T);
b     = 4;                                    % blanket state
h     = 1:3;                                  % external (hidden) states
m     = 5:6;                                  % internal states
N     = 32;                                   % number of probability bins

% get marginals over external and internal states
%--------------------------------------------------------------------------
x{1}  = linspace(-32,32,N);
x{2}  = linspace(-32,32,N);
x{3}  = linspace( 16,64,N);
x{4}  = linspace(-32,32,N);
x{5}  = linspace(-32,32,N);
x{6}  = linspace( 16,64,N);
for i = 1:T
    

    % get marginals over external and internal states
    %----------------------------------------------------------------------
    s    = x;
    s{4} = t(b,i);
    s{5} = 0;
    s{6} = 0;
    
    pH   = spm_softmax(spm_polymtx(s,M.K)*NESS.Ep.Sp);
    pH   = reshape(pH,[N N N]);
    
    s    = x;
    s{1} = 0;
    s{2} = 0;
    s{3} = 0;
    s{4} = t(b,i);
    pI   = spm_softmax(spm_polymtx(s,M.K)*NESS.Ep.Sp);
    pI   = reshape(pI,[N N]);
    
    % conditional (marginal) probabilities
    %----------------------------------------------------------------------
    q     = spm_marginal(pH);
    for j = 1:numel(q)
        
        % marginal
        %------------------------------------------------------------------
        mH(i,:,j) = q{j};
        
        % manifold (conditional expectations)
        %------------------------------------------------------------------
        eH(i,j)   = x{j}*q{j};
    end
    
    q     = spm_marginal(pI);
    for j = 1:numel(q)
        
        % marginal
        %------------------------------------------------------------------
        mI(i,:,j) = q{j};
        
        % manifold (conditional expectations)
        %------------------------------------------------------------------
        eI(i,j)   = x{j}*q{j};
    
    end
    
    q     = spm_marginal(pI);
    for j = 1:numel(q)
        mI(i,:,j) = q{j};
    end
    
    
    eB    = t(b,i);
    
end


%% Notes (symbolic maths or Lorenz system)
%==========================================================================
syms x1 x2 x3 x4 x5 x6 P1 P2 P3 P4 f J x P

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
P    = [P1;P2;P3;P4];
ff   = symfun([P(1)*x(2) - P(1)*(x(1)*(1 - P(4)) + x(4)*P(4));
               P(3)*x(1) - x(2) - x(1)*x(3);
               P(2)*x(3) + x(1)*x(2);
               P(1)*x(5) - P(1)*(x(4)*(1 - P(4)) + x(1)*P(4));
               P(3)*x(4) - x(5) - x(4)*x(6);
               P(2)*x(6) + x(4)*x(5)]/64,x);


J   = [diff(ff,x1) diff(ff,x2) diff(ff,x3) diff(ff,x4) diff(ff,x5) diff(ff,x6)]


%% Notes (generalised coordinates of motion)
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


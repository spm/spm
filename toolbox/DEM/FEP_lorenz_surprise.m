function FEP_lorenz_surprise
% This demo provides an elementary characterisation of stochastic chaos
% using the Lorenz system. Effectively, it uses iterated least-squares to
% solve the Helmholtz decomposition of nonequilibrium steady-state flow
% (i.e., the solution to the Fokker Planck equation) using the Lorentz
% system as an example. This furnishes a generative model for stochastic
% chaos in terms of the underlying potential (log nonequilibrium
% steady-state density) and flow operator, with symmetric and antisymmetric
% (skew symmetric) components. The latter (solenoidal) part of the flow
% operator breaks detailed balance and renders the solution a
% nonequilibrium steady-state (NESS) density.
%
% In virtue of using a polynomial expansion for the nonequilibrium
% potential (i.e., surprisal or self information) one can approximate the
% expected flow with a second order polynomial. This can be regarded as a
% Laplace approximation to the nonequilibrium steady-state density. Further
% constraints can be used to specify the stochastic chaos as (state
% dependent) solenoidal flow around a multivariate Gaussian, which might be
% a reasonable approximation in the setting of random fluctuations.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: FEP_lorenz_surprise.m 8077 2021-03-07 15:44:38Z karl $


%% dynamics and parameters of a Lorentz system (with Jacobian)
%==========================================================================
% dxdt = f(x) + w:  see notes at the end of this script
%--------------------------------------------------------------------------
f    = @(x,v,P,M) [P(1)*x(2) - P(1)*x(1);
                   P(3)*x(1) - x(2) - x(1)*x(3);
                   P(2)*x(3) + x(1)*x(2)]/64;
J    = @(x,v,P,M) [[     -P(1),P(1),     0];
                   [P(3) - x(3), -1, -x(1)];
                   [     x(2), x(1),  P(2)]]/64;
P    = [10; -8/3; 32];                 % parameters
x0   = [1; 1; 24];                     % initial state
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

% solution of stochastic and underlying ordinary differential equation
% -------------------------------------------------------------------------
spm_figure('GetWin','NESS (0)'); clf
T    = 2^10;
U.u  = sparse(T,M(1).m);
U.dt = 1;
t    = spm_int_ode(M.pE,M,U);
r    = spm_int_sde(M.pE,M,U);

subplot(3,1,1)
plot(t), hold on, set(gca,'ColorOrderIndex',1)
plot(r,':'), hold off
title('trajectory','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

subplot(3,1,2)
plot3(t(:,1),t(:,2),t(:,3)), hold on, set(gca,'ColorOrderIndex',1)
plot3(r(:,1),r(:,2),r(:,3),':'), hold off
title('State-space','Fontsize',16)
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
axis square, box off


%% state-space four (Laplace) solution 
%--------------------------------------------------------------------------
N    = 4;                                     % number of bins
d    = 15;
x{1} = linspace(-d,d,N);
x{2} = linspace(-d,d,N);
x{3} = linspace(28 - d,28 + d,N);

% Fokker-Planck operator and equilibrium density
%==========================================================================
NESS = spm_ness_hd(M,x);                      % NESS density
n    = numel(x);                             % dimensionality

% solve for a trajectory using a polynomial approximation to flow 
%--------------------------------------------------------------------------
T     = 1024;
t     = zeros(n,T);
dt    = 1/2;
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
title('trajectory','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

subplot(3,2,4)
plot3(t(1,:),t(2,:),t(3,:))
title('State-space','Fontsize',16)
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
axis square, box off

% marginal nonequilibrium steady-state
%--------------------------------------------------------------------------
N    = 64;
x{1} = linspace(-32,32,N);
x{2} = linspace(-32,32,N);
x{3} = linspace(  0,64,N);
p0   = spm_softmax(spm_polymtx(x)*NESS.Ep.Sp);
p0   = reshape(p0,[N N N]);
X    = spm_ndgrid(x);

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,1)
imagesc(x{3},x{2},squeeze(sum(p0,1)))
hold on, plot(t(3,:),t(2,:),'r.'), hold off
xlabel('3rd state'), ylabel('2bd state')
title('NESS density','Fontsize',14)
axis square xy

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,2)
imagesc(x{3},x{1},squeeze(sum(p0,2)))
hold on, plot(t(3,:),t(1,:),'r.'), hold off
xlabel('3rd state'), ylabel('1st state')
title('NESS density','Fontsize',14)
axis square xy

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,3)
imagesc(x{2},x{1},squeeze(sum(p0,3)))
hold on, plot(t(2,:),t(1,:),'r.'), hold off
xlabel('2nd state'), ylabel('1st state')
title('NESS density','Fontsize',14)
axis square xy

subplot(3,1,3)
for i = 1:n
    plot(NESS.f(:,i),NESS.F(:,i),'.','MarkerSize',1), hold on
end, hold off
title('Flows','Fontsize',16)
xlabel('approximate flow'), ylabel('true flow')
axis square xy, box off


%% illustrate conditional dependencies
%==========================================================================
spm_figure('GetWin','NESS (2)'); clf

% covariance based upon a Laplace approximation to surprisal
%--------------------------------------------------------------------------
% m     = NESS.q0(:)'*NESS.X;
% C     = full(NESS.X'*bsxfun(@times,NESS.q0(:),NESS.X) - m'*m);
[m,C] = spm_ness_cond(n,3,NESS.Ep.Sp);


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
s     = round(linspace(-24,24,5));
for i = 1:numel(s)
    subplot(6,numel(s),numel(s)*2 + i)
    [m,j] = min(abs(x{1} - s(i)));
    imagesc(x{3},x{2},squeeze(p0(j,:,:)).^2)
    axis square xy, axis off
    title('2nd and 3rd','Fontsize',12)
end
s     = round(linspace(-24,24,5));
for i = 1:numel(s)
    subplot(6,numel(s),numel(s)*3 + i)
    [m,j] = min(abs(x{2} - s(i)));
    imagesc(x{3},x{1},squeeze(p0(:,j,:)).^2)
    axis square xy, axis off
    title('1st and 3rd','Fontsize',12)
end
s     = round(linspace(16,48,5));
for i = 1:numel(s)
    subplot(6,numel(s),numel(s)*4 + i)
    [m,j] = min(abs(x{2} - s(i)));
    imagesc(x{2},x{1},squeeze(p0(:,:,j)).^2)
    axis square xy, axis off
    title('1st and 2nd','Fontsize',12)
end


%% illustrate conditional distributions
%==========================================================================

% get marginals over external and internal states
%--------------------------------------------------------------------------
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
   
    % manifold (conditional expectations)
    %-------------------------------------------====-----------------------
    eH(i)   = x{h}*pH(:,i);
    
end

subplot(6,1,6), hold off
imagesc(1:T,x{h},(1 - pH).^2), hold on
plot(t(h,:),'r'), plot(eH,'w'), axis xy
title('Conditional density and expectations','Fontsize',14)
xlabel('time'), ylabel('1st state')


%% functional form of polynomial flow
%==========================================================================:
clear Ep
syms  W 'real'
syms  x [1 n] 'real'
sympref('FloatingPointOutput',0);

% symbolic coefficients of potential
%--------------------------------------------------------------------------
o     = NESS.o;
nb    = size(NESS.o,2);
for i = 1:nb
    o     = NESS.o(:,i);
    if sum(o) == 0
        str = 's0';
    elseif sum(o) == 1
        str = sprintf('s%i',find(o));
    else
        k   = find(o);
        if numel(k) < 2
        str = sprintf('s%i%i',k,k);
        else
            str = sprintf('s%i%i',k(1),k(2));
        end
    end
    eval([str '= NESS.Ep.Sp(i);']);
    Ep.Sp(i,1) = sym(str,'real');
    if abs(NESS.Ep.Sp(i)) < 1e-4
        Ep.Sp(i,1) = 0;
    end

end

% symbolic coefficients of flow operator
%--------------------------------------------------------------------------
k     = 0;
for j = 1:(n*n - n)/2
    for i = 1:nb
        k = k + 1;
        o     = NESS.o(:,i);
        if sum(o) == 0
            str = 'q0';
        elseif sum(o) == 1
            str = sprintf('q%i',find(o));
        else
            d   = find(o);
            if numel(d) < 2
                str = sprintf('q%i%i',d,d);
            else
                str = sprintf('q%i%i',d(1),d(2));
            end
        end
        eval([str '= NESS.Ep.Qp(k);']);
        Ep.Qp(k,1) = sym(str,'real');
        if abs(NESS.Ep.Qp(k)) < 1e-4
            Ep.Qp(k,1) = 0;
        end
    end
end


% replace sample points with symbolic variables
%--------------------------------------------------------------------------
M.X  = x;
M.W  = diag(kron(ones(n,1),W));

% evaluate flow, flow operators and Hessians for display
%--------------------------------------------------------------------------
U           = spm_ness_U(M);
[F,Q,S,L,H] = spm_NESS_gen(Ep,M);
F = F'
Q = reshape(cat(1,Q{:}),n,n)
S = S
L = L'
H = H
f = real(U.f)

disp('f = ')
disp(latex(f)), disp(' ')
disp('f = ')
disp(latex(F)), disp(' ')
disp('Q = ')
disp(latex(Q)), disp(' ')
disp('L = ')
disp(latex(L)), disp(' ')
disp('H = ')
disp(latex(H)), disp(' ')

% ensure the Laplace approximation is chaotic
%==========================================================================

% Jacobian function with polynomial form
%--------------------------------------------------------------------------
J    = matlabFunction([diff(F,x1) diff(F,x2) diff(F,x3)]);

% state space for evaluation
%--------------------------------------------------------------------------
N    = 16;
d    = 23;
x    = cell(3,1);
x{1} = linspace(-d,d,N);
x{2} = linspace(-d,d,N);
x{3} = linspace(28 - d,28 + d,N);
p0   = spm_softmax(spm_polymtx(x)*NESS.Ep.Sp);
p0   = reshape(p0,[N N N]);
X    = spm_ndgrid(x);
W    = 1/16;

% evaluate eigenvalues of Jacobian
%--------------------------------------------------------------------------
nX    = size(X,1);
E     = zeros(n,nX);
for i = 1:nX
    s      = full(X(i,:));
    E(:,i) = eig(J(W,q0,q1,q2,q3,q11,q12,q13,q22,q23,q33,s3,s11,s12,s22,s33,s(1),s(2),s(3)));
end

% Lyapunov exponent and Hausdorff dimension (Kaplan-Yorke conjecture)
%--------------------------------------------------------------------------
LE  = sort(real(spm_dot(E,p0(:))),'descend')
j   = sum(LE > 0);
CD  = j + sum(LE(1:j))/abs(LE(j + 1))

clear syms




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

P    = [10; -8/3; 32; -1/2];            % parameters
x0   = [1; 1; 32; 1; 1; 30];            % initial states
W    = diag([1 1 1 1 1 1]/16);          % precision of random fluctuations

% state space model
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;
M.W  = W;

% state space
%--------------------------------------------------------------------------
N    = 3;                                     % number of bins
d    = 12;
x{1} = linspace(-d,d,N);
x{2} = linspace(-d,d,N);
x{3} = linspace(28 - d,28 + d,N);
x{4} = linspace(-d,d,N);
x{5} = linspace(-d,d,N);
x{6} = linspace(28 - d,28 + d,N);

% Fokker-Planck operator and equilibrium density
%==========================================================================
NESS = spm_ness_hd(M,x);                      % NESS density
n    = numel(x);                              % dimensionality

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

% illustrate convergence to nonequilibrium steady-state
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (3)'); clf
subplot(3,2,3)
plot(t')
title('Trajectory','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

subplot(3,2,4)
plot3(t(2,:),t(1,:),t(3,:))
hold on, plot3(t(5,:),t(4,:),t(6,:)), hold off
title('State-space','Fontsize',16)
xlabel('2nd states'), ylabel('1st states'), zlabel('3rd states')
axis square, box off

% marginal nonequilibrium steady-state
%--------------------------------------------------------------------------
N    = 32;
s    = cell(6,1);
s{1} = linspace(-32,32,N);
s{2} = linspace(-32,32,N);
s{3} = linspace(  0,64,N);
s{4} = 0;
s{5} = 0;
s{6} = 0;
p0   = spm_softmax(spm_polymtx(s)*NESS.Ep.Sp);
p0   = reshape(p0,[N N N]);

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,1)
imagesc(s{3},s{2},squeeze(sum(p0,1)))
hold on, plot(t(3,:),t(2,:),'r.'), hold off
xlabel('3rd state'), ylabel('2nd state')
title('NESS density','Fontsize',14)
axis square xy

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,2)
imagesc(s{3},s{1},squeeze(sum(p0,2)))
hold on, plot(t(3,:),t(1,:),'r.'), hold off
xlabel('3rd state'), ylabel('1st state')
title('NESS density','Fontsize',14)
axis square xy

% nonequilibrium steady-state (marginal)
%--------------------------------------------------------------------------
subplot(3,3,3)
imagesc(s{2},s{1},squeeze(sum(p0,3)))
hold on, plot(t(2,:),t(1,:),'r.'), hold off
xlabel('2nd state'), ylabel('1st state')
title('NESS density','Fontsize',14)
axis square xy

subplot(3,1,3)
for i = 1:n
    plot(NESS.f(:,i),NESS.F(:,i),'.','MarkerSize',1), hold on
end, hold off
title('Flows','Fontsize',16)
xlabel('approximate flow'), ylabel('true flow')
axis square xy, box off


%% illustrate conditional dependencies
%==========================================================================
spm_figure('GetWin','NESS (4)'); clf

% mean and covariance
%--------------------------------------------------------------------------
[m,C] = spm_ness_cond(n,3,NESS.Ep.Sp);

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


%% illustrate extrinsic information geometry
%==========================================================================

% get conditional marginals over external and internal states
%--------------------------------------------------------------------------
b     = 1;                                     % blanket state
E     = zeros(n - 1,T);
V     = zeros(n - 1,T);
for i = 1:T
    [m,C]  = spm_ness_cond(n,3,NESS.Ep.Sp,b,t(b,i));
    E(:,i) = m;
    V(:,i) = diag(C);
end

subplot(4,1,2)
plot(t'), title('Coupled dynamics','Fontsize',14)
ylabel('states'), spm_axis tight, box off

subplot(4,1,3)
spm_plot_ci(E(4,:),V(4,:)')
title('Expected internal state','Fontsize',14)
ylabel('states'), spm_axis tight, box off

subplot(4,1,4), hold off
spm_plot_ci(E(1,:),V(1,:)'), hold on
plot(t(2,:),'-.'), hold off
title('Expected external state','Fontsize',14)
ylabel('states'), spm_axis tight, box off


%% Notes for functional form of polynomial flow
%--------------------------------------------------------------------------
clear Ep
syms  W  'real'
syms  x  [1 n] 'real'
syms  bs 'real'

% set small polynomial coefficients to zero
%--------------------------------------------------------------------------
o     = NESS.o;
for i = 1:numel(NESS.Ep.Sp)
    o     = NESS.o(:,i);
    if sum(o) == 0
        str = 's0';
    elseif sum(o) == 1
        str = sprintf('s%i',find(o));
    else
        k   = find(o);
        if numel(k) < 2
        str = sprintf('s%i%i',k,k);
        else
            str = sprintf('s%i%i',k(1),k(2));
        end
    end
    Ep.Sp(i,1) = sym(str,'real');
    if abs(NESS.Ep.Sp(i)) < 1e-4
        Ep.Sp(i,1) = 0;
    end
end

% display conditional moments
%--------------------------------------------------------------------------
sympref('FloatingPointOutput',0);
[m,C]  = spm_ness_cond(n,3,Ep.Sp,b,bs)

disp('/nmu = ')
disp(latex(m(4))),      disp(' ')
disp('E[q] = ')
disp(latex(m(1))),      disp(' ')
disp('Cov[q] = ')
disp(latex(C(1,1))),    disp(' ')
disp('E[h] = E[m] x ')
disp(latex(m(1)/m(4))), disp(' ')


return


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


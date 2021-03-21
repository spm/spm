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
% $Id: FEP_lorenz_surprise.m 8085 2021-03-21 12:27:26Z karl $


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
P    = [10; -8/3; 32];                 % parameters [sig, beta, rho]
x0   = [1; 1; 24];                     % initial state
W    = diag([1 1 1]/16);               % precision of random fluctuations

disp('Lyapunov dimnesion')
disp(3 - 2*(P(1) - P(2) + 1)/(P(1) + 1 + sqrt((P(1) - 1)^2 + 4*P(1)*P(3))))

% state-space model (for SPM integrators)
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;
M.W  = W;
M.K  = 3;

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
title('Trajectory','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

subplot(3,1,2)
plot3(t(:,1),t(:,2),t(:,3)), hold on, set(gca,'ColorOrderIndex',1)
plot3(r(:,1),r(:,2),r(:,3),':'), hold off
title('State-space','Fontsize',16)
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
axis square, box off

% priors over potential parameters
%--------------------------------------------------------------------------
% M.pS = spm_ness_m2S(mean(r),cov(r));


%% state-space four (Laplace) solution 
%--------------------------------------------------------------------------
N    = 4;                                     % number of bins
d    = 18;
m    = [0 0 28];
x{1} = linspace(m(1) - d,m(1) + d,N);
x{2} = linspace(m(2) - d,m(2) + d,N);
x{3} = linspace(m(3) - d,m(3) + d,N);

% Fokker-Planck operator and equilibrium density
%==========================================================================
NESS = spm_ness_hd(M,x);                       % NESS density
n    = numel(x);                               % dimensionality

%% functional form of polynomial flow
%==========================================================================
Ep    = NESS.Ep;
syms  w 'real'
syms  x [1 n] 'real'
syms  s [numel(Ep.Sp) 1] 'real'
syms  q [numel(Ep.Qp) 1] 'real'

% equation of motion as a Matlab function
%--------------------------------------------------------------------------
Ep.Sp   = s;
Ep.Qp   = q;
O       = M;
O.X     = x;
O.W     = diag(kron(ones(n,1),w));
[fs,ss,qs,ls,hs] = spm_NESS_gen(Ep,O);

qs      = reshape(cat(1,qs{:}),n,n);
M.fs    = matlabFunction(fs);
M.ss    = matlabFunction(ss);
M.qs    = matlabFunction(qs);
M.ls    = matlabFunction(ls);
M.hs    = matlabFunction(hs);


%% evaluate flow as a Matlab function of states
%--------------------------------------------------------------------------
w   = {W(1)};
Sp  = num2cell(NESS.Ep.Sp);
Qp  = num2cell(NESS.Ep.Qp);

sympref('FloatingPointOutput',1);
disp('ss = ')
disp(latex(M.ss(Sp{:},x(1),x(2),x(3))')), disp(' ')
disp('qs = ')
disp(latex(M.qs(Qp{:},w{:},x(1),x(2),x(3))')), disp(' ')
disp('fs = ')
disp(latex(M.fs(Qp{:},Sp{2:end},w{:},x(1),x(2),x(3))')), disp(' ')

% solve for a trajectory using a polynomial approximation to flow 
%--------------------------------------------------------------------------
T     = 4096;
t     = zeros(n,T);
dt    = 1/4;
s     = x0;
for i = 1:T
    x      = num2cell(s');
    ds     = M.fs(Qp{:},Sp{2:end},w{:},x{:});
    s      = s + ds'*dt; 
    t(:,i) = s;                                % expected trajectory
end

% illustrate convergence to nonequilibrium steady-state
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (1)'); clf
subplot(3,1,1)
plot(t')
title('trajectory','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off


%% solve for a trajectory using a polynomial approximation to flow 
%--------------------------------------------------------------------------
T     = 2048;
t     = zeros(n,T);
u     = zeros(n,T);
dt    = 1/4;
s     = x0;
r     = x0;
Ep    = NESS.Ep;
for i = 1:T
    M.X    = s';
    [ds,j] = spm_NESS_gen(Ep,M);
    s      = s + ds'*dt; 
    t(:,i) = s;                                % expected trajectory
    LF(i)  = j;                                % Lyapunvov function
    
    M.X    = r';
    dr     = spm_NESS_gen(Ep,M);
    r      = r + dr'*dt + sqrtm(inv(W))*randn(n,1)/2*dt; 
    u(:,i) = r;                                % stochastic solution
end

% illustrate convergence to nonequilibrium steady-state
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (1)'); clf
subplot(3,2,3)
plot(t'),     hold on, set(gca,'ColorOrderIndex',1)
plot(u',':'), hold off
title('trajectory','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

subplot(3,2,4)
plot3(t(1,:),t(2,:),t(3,:)),     hold on, set(gca,'ColorOrderIndex',1)
plot3(u(1,:),u(2,:),u(3,:),':'), hold on
title('State-space','Fontsize',16)
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
axis square, box off

% marginal nonequilibrium steady-state
%--------------------------------------------------------------------------
N    = 64;
x{1} = linspace(-32,32,N);
x{2} = linspace(-32,32,N);
x{3} = linspace(  0,64,N);
p0   = spm_softmax(spm_polymtx(x,M.K)*NESS.Ep.Sp);
p0   = reshape(p0,[N N N]);

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

subplot(3,2,5)
plot(LF)
title('Potential','Fontsize',16)
xlabel('time'), ylabel('self-information')
axis square xy, box off

subplot(3,2,6)
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
% m   = NESS.q0(:)'*NESS.X;
% C   = full(NESS.X'*bsxfun(@times,NESS.q0(:),NESS.X) - m'*m);
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
    %----------------------------------------------------------------------
    eH(i)   = x{h}*pH(:,i);
    
end

subplot(6,1,6), hold off
imagesc(1:T,x{h},(1 - pH).^2), hold on
plot(t(h,:),'r'), plot(eH,'w'), axis xy
title('Conditional density and expectations','Fontsize',14)
xlabel('time'), ylabel('1st state')


%% functional form of polynomial flow
%==========================================================================
TOL   = 1e-4;
Ep    = NESS.Ep;
syms  w 'real'
syms  x [1 n] 'real'
syms  s [numel(Ep.Sp) 1] 'real'
syms  q [numel(Ep.Qp) 1] 'real'
sympref('FloatingPointOutput',0);

Sp  = NESS.Ep.Sp;
Qp  = NESS.Ep.Qp;
i   = abs(Sp) < TOL; s(i) = 0;
i   = abs(Qp) < TOL; q(i) = 0;

% replace sample points with symbolic variables
%--------------------------------------------------------------------------
Ep.Sp = s;
Ep.Qp = q;
O     = M;
O.X   = x;
O.W   = diag(kron(ones(n,1),w));

% evaluate flow, flow operators and Hessians for display
%--------------------------------------------------------------------------
U           = spm_ness_U(O);
[F,S,Q,L,H] = spm_NESS_gen(Ep,O);
F = F'
Q = reshape(cat(1,Q{:}),n,n)
S = S
L = L'
H = H
D = [diff(S,x1); diff(S,x2); diff(S,x3)];
f = real(U.f);

disp('f = ')
disp(latex(f)), disp(' ')
disp('F = ')
disp(latex(F)), disp(' ')
disp('Q = ')
disp(latex(Q)), disp(' ')
disp('L = ')
disp(latex(L)), disp(' ')
disp('H = ')
disp(latex(H)), disp(' ')

disp('S = ')
disp(latex(S)), disp(' ')
disp('dS/dx = ')
disp(latex(D)), disp(' ')

% ensure the Laplace approximation is chaotic
%==========================================================================

% Jacobian function with polynomial form
%--------------------------------------------------------------------------
syms  w 'real'
syms  x [1 n] 'real'
syms  s [numel(Ep.Sp) 1] 'real'
syms  q [numel(Ep.Qp) 1] 'real'
sympref('FloatingPointOutput',1);

J    = symfun([diff(F,x1), diff(F,x2), diff(F,x3)],[q',s',w,x]);
J    = matlabFunction(J);

% state space for evaluation
%--------------------------------------------------------------------------
N    = 28;
d    = 28;
x    = cell(3,1);
x{1} = linspace(-d,d,N);
x{2} = linspace(-d,d,N);
x{3} = linspace(28 - d,28 + d,N);
p0   = spm_softmax(spm_polymtx(x)*NESS.Ep.Sp);
p0   = reshape(p0,[N N N]);
X    = spm_ndgrid(x);
w    = W(1);

% evaluate eigenvalues of Jacobian
%--------------------------------------------------------------------------
Sp  = num2cell(NESS.Ep.Sp);
Qp  = num2cell(NESS.Ep.Qp);
nX  = size(X,1);
E   = zeros(n,nX);
Ji  = zeros(n,n,nX);
for i = 1:nX
    s         = num2cell(full(X(i,:)));
    Ji(:,:,i) = J(Qp{:},Sp{:},w,s{:});
    E(:,i)    = sort(real(eig(Ji(:,:,i))),'descend');
end

% Lyapunov exponent and Hausdorff dimension (Kaplan-Yorke conjecture)
%--------------------------------------------------------------------------
LE  = spm_dot(E,p0(:))
j   = sum(LE > 0);
CD  = j + sum(LE(1:j))/abs(LE(j + 1))

clear


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
x0   = [1; 1; 32; 1; 0; 24];            % initial states
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
N    = 4;                                     % number of bins
d    = 18;
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

%% solve for a trajectory 
%--------------------------------------------------------------------------
T     = 1024;
t     = zeros(n,T);
u     = zeros(n,T);
dt    = 1/2;
s     = x0;
r     = x0;
Ep    = NESS.Ep;
for i = 1:T
    M.X    = s';
    [ds,j] = spm_NESS_gen(Ep,M);
    s      = s + ds'*dt; 
    t(:,i) = s;                                % expected trajectory
    LF(i)  = j;                                % Lyapunvov function
    
    M.X    = r';
    dr     = spm_NESS_gen(Ep,M);
    r      = r + dr'*dt + sqrtm(inv(W))*randn(n,1)/2*dt; 
    u(:,i) = r;                                % stochastic solution
end

% illustrate convergence to nonequilibrium steady-state
%--------------------------------------------------------------------------
spm_figure('GetWin','NESS (3)'); clf
subplot(3,2,3)
plot(t'),     hold on, set(gca,'ColorOrderIndex',1)
plot(u',':'), hold off
title('Trajectory','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

subplot(3,2,4), hold off
plot3(t(5,:),t(4,:),t(6,:),'c'), hold on
plot3(t(2,:),t(1,:),t(3,:),'b')
plot3(u(5,:),u(4,:),u(6,:),':c')
plot3(u(2,:),u(1,:),u(3,:),':b')

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
TOL   = 1e-4;
Sp    = NESS.Ep.Sp;
i     = abs(Sp) < TOL; Sp(i) = 0;
[m,C] = spm_ness_cond(n,3,Sp);

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
plot(E(4,:))
title('Conditional expectation of internal state','Fontsize',14)
ylabel('state'), spm_axis tight, box off

subplot(4,1,4), hold off
spm_plot_ci(E(1,:),V(1,:)'), hold on
plot(t(2,:),'-.'), hold off
title('Conditional density over external state','Fontsize',14)
ylabel('states'), spm_axis tight, box off


%% functional form of conditional independencies
%--------------------------------------------------------------------------
clear Ep
syms  W  'real'
syms  x  [1 n] 'real'
syms  bs 'real'

% set small polynomial coefficients to zero
%--------------------------------------------------------------------------
TOL   = 1e-4;
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
    if abs(NESS.Ep.Sp(i)) < TOL
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

%% functional form of polynomial flow
%==========================================================================
Ep    = NESS.Ep;
syms  w 'real'
syms  x [1 n] 'real'
syms  s [numel(Ep.Sp) 1] 'real'
syms  q [numel(Ep.Qp) 1] 'real'
sympref('FloatingPointOutput',0);

TOL   = 1e-4;
i     = abs(NESS.Ep.Sp) < TOL; s(i) = 0;
i     = abs(NESS.Ep.Qp) < TOL; q(i) = 0;

% replace sample points with symbolic variables
%--------------------------------------------------------------------------
Ep.Sp = s;
Ep.Qp = q;
O     = M;
O.X   = x;
O.W   = diag(kron(ones(n,1),w));

% evaluate flow, flow operators and Hessians for display
%--------------------------------------------------------------------------
[F,S,Q,L,H] = spm_NESS_gen(Ep,O);

for i = 1:n
    for j = i:n
        fprintf('Q = %i%i\n',i,j)
        disp(Q{1,1}), disp(latex(Q{1,1})), disp(' ')
    end
end
disp('L = ')
disp(latex(L')), disp(' ')
disp('H = ')
disp(latex(H)), disp(' ')

return


function varargout = num2csl(a)
varargout = num2cell(a);


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


f
disp(latex(f)),  disp(' ')
ff
disp(latex(ff)), disp(' ')


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


%% notes to identify parameters for flow constraints
%--------------------------------------------------------------------------
[b,D,H,o] = spm_polymtx([{0},{0},{0}],3);
n    = numel(D);                               % dimensionality
nb   = numel(b);
DDG  = 0;
if DDG
    nu = (n^2 - n)/2;
else
    nu = (n^2 + n)/2;
end
nB   = nu*nb;

% dxdt = f(x) + w:  see notes at the end of this script
%--------------------------------------------------------------------------
f    = @(x,v,P,M) [P(1)*x(2) - P(1)*x(1);
                   P(3)*x(1) - x(2) - x(1)*x(3);
                   P(2)*x(3) + x(1)*x(2)]/64;
J    = @(x,v,P,M) [[     -P(1),P(1),     0];
                   [P(3) - x(3), -1, -x(1)];
                   [     x(2), x(1),  P(2)]]/64;
P    = [10; -8/3; 32];                 % parameters [sig, beta, rho]
x0   = [1; 1; 24];                     % initial state


% state-space model (for SPM integrators)
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;


% constraints on flow operator (given definiteness of Hessian)
%-------------------------------------------------------------------------- 
syms  w 'real'
syms  x [1 n] 'real'
syms  s [nb 1] 'real'
syms  q [nB 1] 'real'
M.X   = x;
M.W   = diag(kron(ones(n,1),w));

i     = find(~any(U.o == 2)); s(i) = 0;
J     = M.J(ones(1,n),[],P,M);
Ep.Sp = s;
Ep.Qp = q;
F     = spm_NESS_gen(Ep,M);
iQ    = [];
for i = 1:n
    for j = 1:n
        if ~J(i,j)
            Js(i,j)   = diff(F(i),x(j));
            for k = 1:numel(q)
                kQ(k) = diff(Js(i,j),q(k));
            end
            iQ   = unique([iQ,find(kQ)]);
        end
    end
end

% constraints on flow operator (given definiteness of G)
%--------------------------------------------------------------------------
syms  s [nb 1] 'real'
syms  q [nB 1] 'real'
Ep.Sp = s;
Ep.Qp = q*0;
iS    = [];
F     = spm_NESS_gen(Ep,M);
for i = 1:n
    for j = 1:n
        if ~J(i,j)
            Jq(i,j)   = diff(F(i),x(j));
            for k = 1:numel(s)
                kS(k) = diff(Jq(i,j),s(k));
            end
            iS   = unique([iS,find(kS)]);
        end
    end
end


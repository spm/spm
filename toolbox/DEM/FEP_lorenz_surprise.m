function FEP_lorenz_surprise
% This demo computes the cost-function (negative reward) for a Lorenz
% system; to show cost can be computed easily from value (negative 
% surprise or sojourn time). However, value is not a Lyapunov function
% because the flow is not curl-free (i.e., is not irrotational).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: FEP_lorenz_surprise.m 8000 2020-11-03 19:04:17Z karl $


% generative model
%==========================================================================                       % switch for demo
% f   = @(x,v,P,M) [-P(1) 0 0; 0 -P(2) 0; 0 0 -P(3)]*x/64;
% P   = [1; 1; 8];
% x0  = [0; 0; 0];
% W   = diag([1; 1; 1/4]);
 
% dynamics and parameters
%--------------------------------------------------------------------------
f    = @(x,v,P,M) [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64;
J    = @(x,v,P,M) [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]/64 + ...
                  [0 0 0; -x(3) 0 0; 0 x(1) 0]/64;
P    = [10; -8/3; 32];
x0   = [1; 1; 32];
W    = diag([1 1 1]/16);

% state space model
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;
M.W  = W;                           % precision of State noise

% state space
%--------------------------------------------------------------------------
N    = 16;
x{1} = linspace(-32,32,N);
x{2} = linspace(-32,32,N);
x{3} = linspace( 16,48,N);

% Fokker-Planck operator and equilibrium density
%==========================================================================
[q0,X,F,f,H,J] = spm_ness_hd(M,x);
nx             = size(q0);
n              = numel(nx);

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


% illustrate convergence to nonequilibrium steady-state
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM'); clf
subplot(3,2,5)
plot(t')
axis square
title('trajectory','Fontsize',16)

subplot(3,2,6)
plot3(t(1,:),t(2,:),t(3,:))
axis square
title('trajectory','Fontsize',16)

% axes and trajectory
% --------------------------------------------------------------------------
% T    = 2^13;
% U.u  = sparse(T,M(1).m);
% U.dt = 1;
% t    = spm_int_sde(M.pE,M,U);




% nonequilibrium steady-state
%--------------------------------------------------------------------------
subplot(3,2,1)
imagesc(x{3},x{2},exp(log(squeeze(sum(q0,1)))/8))
hold on, plot(t(3,:),t(2,:),'r'), hold off
axis square xy
title('value','Fontsize',16)

% nonequilibrium steady-state
%--------------------------------------------------------------------------
subplot(3,2,2)
imagesc(x{3},x{1},exp(log(squeeze(sum(q0,2)))/8))
hold on, plot(t(3,:),t(1,:),'r'), hold off
axis square xy
title('value','Fontsize',16)

% nonequilibrium steady-state
%--------------------------------------------------------------------------
subplot(3,2,3)
imagesc(x{2},x{1},exp(log(squeeze(sum(q0,3)))/8))
hold on, plot(t(2,:),t(1,:),'r'), hold off
axis square xy
title('value','Fontsize',16)

subplot(3,2,4)
for i = 1:3
    plot(f(i,:)',F(:,i),'.','MarkerSize',1), hold on
end, hold off
axis square xy
title('value','Fontsize',16)


return

% nonequilibrium steady-state
%--------------------------------------------------------------------------
subplot(3,2,1)
imagesc(x{3},x{2},squeeze(sum(q0,1)))
hold on, plot(t(:,3),t(:,2),'r'), hold off
axis square xy
title('value','Fontsize',16)

% nonequilibrium steady-state
%--------------------------------------------------------------------------
[d,z] = min(abs(x{1} - M.x(1)));
subplot(3,2,2)
m = squeeze(q0(z,:,:));
imagesc(x{3},x{2},m)
hold on
i = find(abs(t(:,1) - x{1}(z)) < 4);
plot(t(:,3),t(:,2),'r')
plot(t(i,3),t(i,2),'.c','MarkerSize',1)
hold off
axis square xy
title('value','Fontsize',16)


% gradients based upon Helmholtz decomposition
%==========================================================================

% state space
%--------------------------------------------------------------------------
N    = 32;
x{1} = linspace(-32,32,N);
x{2} = linspace(-32,32,N);
x{3} = linspace(0,64,N);

[q0,X,H,J] = spm_ness_hd(M,x);


% covariance
%--------------------------------------------------------------------------
% m    = q0(:)'*X;
% full(X'*bsxfun(@times,q0(:),X) - m'*m);

% Gradients of surprisal at nonequilibrium steady-state
%--------------------------------------------------------------------------
subplot(3,2,4)
imagesc(x{3},x{2},exp(log(squeeze(sum(q0,1)))/2))
axis square xy
title('cost','Fontsize',16)


% Gradients of surprisal at nonequilibrium steady-state
%--------------------------------------------------------------------------
subplot(3,4,5)
imagesc(J), axis square
title('Jacobian','Fontsize',12)

subplot(3,4,6)
imagesc(H), axis square
title('Hessian','Fontsize',12)


N     = 32;
x{2}  = linspace(-32,32,N);
x{3}  = linspace(0,64,N);
for b = linspace(-20,20,32)
    
    x{1} = linspace(-32,32,9) + b;
    q    = spm_ness_hd(M,x);
    q    = squeeze(q(4,:,:));
    
    subplot(3,2,3)
    imagesc(x{3},x{2},squeeze(exp(log(q)/2)))
    axis square xy
    title('Conditional','Fontsize',12)
    drawnow,disp(b)
    
end












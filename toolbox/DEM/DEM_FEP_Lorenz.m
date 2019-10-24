function DEM_FEP_Lorenz
%--------------------------------------------------------------------------
% This is a simple demonstration of deterministic convergence to
% nonequilibrium steady-state, using the Lorenz system. Deterministic
% solutions (with a Rayleigh parameter of 28) are obtained for 2048 initial
% states, integrating over eight seconds (with a time step of 1/64
% seconds). Crucially, the initial autonomous states are the same for each
% solution and yet the final density over the autonomous (i.e., active)
% state converges to the non-equilibrium steady-state density over time.
% This is apparent in the collapse of the divergence between the sample
% densities (over all states) and the final (NESS) density - as evaluated
% simply using a Gaussian approximation to the ensemble densities at each
% point in time. The upper plots show the propagated states at four points
% in time. As time progresses, this density comes to assume the familiar
% butterfly form of the Lorenz attractor. However, these states are not
% trajectories through state space, they are the endpoints of paths from an
% ensemble of starting locations (shown in the right plot). In this
% illustration, we treat the first state of the Lorenz system as the active
% state, the second state as the sensory state and the third state plays
% the role of an external or hidden state. This designation is based upon
% the fact that the first state is not influenced by the first. In short,
% this numerical example shows how uncertainty about external states is
% propagated over time to induce uncertainty about a particle's state; even
% when the initial (particular) state is known.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_FEP_Lorenz.m 7679 2019-10-24 15:54:07Z spm $

% generative model
%==========================================================================                       % switch for demo
spm_figure('GetWin','DEM'); clf

% flow
%--------------------------------------------------------------------------
G.f   = @(x,v,P,G) v(:) + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64;

% set up
%--------------------------------------------------------------------------
ax    = [-20 20 -20 60];              % axis range
dt    = 1/64;                         % time step
T     = 512;                          % length of trajectory
N     = 2^13;                         % number of paths 
P     = [10; -8/3; 28];               % Rayleigh parameter (28 or 8)
g     = exp(0);                      % random fluctuations

% initial states
%--------------------------------------------------------------------------
U.u   = randn(T,3)*g;
G.x   = [1; 1; 24];
x     = spm_int_L(P,G,U);
i     = T/4:T;
x0    = mean(x(i,:))';
s0    = std(x(i,3));
for k = 1:N
    
    % integrate timeseries with random initial hidden state
    %----------------------------------------------------------------------
    U.u      = randn(T,3)*g;    % random fluctuations
    G.x      = x0 + [0;0;randn*s0];
    x(:,:,k) = spm_int_L(P,G,U);
    
end

% final (NESS) Gaussian density
%--------------------------------------------------------------------------
m     = mean(squeeze(x(T,:,:))');
c     = cov(squeeze(x(T,:,:))');

% convergence to nonequilibrium state
%--------------------------------------------------------------------------
tt    = 1:(T - 1);
for t = 1:numel(tt)
    
    % divergence from nonequilibrium steady-state
    %----------------------------------------------------------------------
    mt   = mean(squeeze(x(tt(t),:,:))');
    ct   = cov(squeeze(x(tt(t),:,:))');
    D(t) = spm_kl_normal(mt,ct,m,c);
    
    
    % information length
    %----------------------------------------------------------------------
    mdt   = mean(squeeze(x(tt(t) + 1,:,:))');
    cdt   = cov(squeeze(x(tt(t) + 1,:,:))');
    dL(t) = spm_kl_normal(mdt,cdt,mt,ct);
    
end

% information length
%--------------------------------------------------------------------------
dL = dL.*(dL > exp(-8));
L  = cumsum(sqrt(dL));
L  = L - L(end);

% plot ensemble densities
%--------------------------------------------------------------------------
td    = fix(linspace(1,T,4));
b1    = linspace(ax(1),ax(2),128);
b2    = linspace(ax(3),ax(4),128);
d1    = b1(2) - b1(1);
d2    = b2(2) - b2(1);
for t = 1:4
    
    % sample density
    %----------------------------------------------------------------------
    xt    = squeeze(x(td(t),:,:))';
    p     = zeros(numel(b1),numel(b2)) + 1;
    for i = 1:N
        i1  = round((xt(i,1) - ax(1))/d1);
        i2  = round((xt(i,3) - ax(3))/d2);
        try
            p(i1,i2) = p(i1,i2) + 1;
        end
    end
    p     = p/sum(p(:));
    
    % plot
    %----------------------------------------------------------------------
    subplot(3,4,t),imagesc(b1,b2,1- spm_conv(p',2,2))
    axis square, axis xy
    title('Predictive density')
    xlabel('active state'),ylabel('external state'),drawnow
end

% and Kulback-Leibler divergence
%--------------------------------------------------------------------------
subplot(3,1,2)
plot(tt*dt,D,tt*dt,L,':b',tt*dt,D*0,'-.b')
title('Kullback-Leibler divergence','Fontsize',16)
xlabel('time (seconds)'),ylabel('divergence (nats)'), box off
set(gca,'YLim',[-1.2,1.2]*min(D(1),128))

legend('KL divergence','Information length'),legend BOXOFF

return

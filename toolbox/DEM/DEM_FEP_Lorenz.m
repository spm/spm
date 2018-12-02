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
% densities (over all states) and the final (NESS) density – as evaluated
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
% propagated over time to induce uncertainty about a particle’s state; even
% when the initial (particular) state is known.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_FEP_Lorenz.m 7502 2018-12-02 12:28:03Z karl $

% generative model
%==========================================================================                       % switch for demo
spm_figure('GetWin','DEM'); clf

% flow
%--------------------------------------------------------------------------
G.f   = @(x,v,P,G) v(:) + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64;

% set up
%--------------------------------------------------------------------------
dt    = 1/64;                       % time step
T     = 512;                        % length of trajectory
N     = 2048;                       % number of paths
U.u   = zeros(T,3);                 % inputs
P     = [10; -8/3; 28];             % Rayleigh parameter
for k = 1:N
    
    % integrate timeseries with random initial hidden state
    %----------------------------------------------------------------------
    G.x      = [1; 1; 25 + randn*4];
    x(:,:,k) = spm_int_L(P,G,U);
    
end

% final (NESS) Gaussian density
%--------------------------------------------------------------------------
m     = mean(squeeze(x(T,:,:))');
c     = cov(squeeze(x(T,:,:))');

% convergence to nonequilibrium state
%--------------------------------------------------------------------------
D     = [];
tt    = fix(linspace(1,T,64));
for t = 1:numel(tt)
    mt   = mean(squeeze(x(tt(t),:,:))');
    ct   = cov(squeeze(x(tt(t),:,:))');
    D(t) = spm_kl_normal(mt,ct,m,c);
end

% plot ensemble densities
%--------------------------------------------------------------------------
td    = fix(linspace(1,T,4));
for t = 1:4
    xt = squeeze(x(td(t),:,:))';
    subplot(3,4,t),plot(xt(:,1),xt(:,3),'.b','MarkerSize',1)
    axis square, axis([-20 20 0 60])
    title('Predictive density')
    xlabel('active state'),ylabel('external state'),drawnow
    
end

% and Kulback-Leibler divergence
%--------------------------------------------------------------------------
subplot(3,1,2)
plot(tt*dt,D)
title('Kulback-Leibler divergence','Fontsize',16)
xlabel('time (seconds)'),ylabel('divergence (nats)')

return

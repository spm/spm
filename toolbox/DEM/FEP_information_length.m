function FEP_information_length(gi,qi,ci)
%__________________________________________________________________________
% This demonstration routine illustrates the key role of solenoidal flow
% (that breaks detailed balance) in optimisation and self-organisation. The
% first section shows that increasing the solenoid flow leads to mixing
% that accelerates the convergence of a random dynamical system to its
% (nonequilibrium) steady-state or minimum. Heuristically, solenoidal flow
% —on the level set of on objective function (here the log density of the
% said steady state)—can be regarded as searching for ‘points of entry’
% in state space with steep gradients. The key observation here is that the
% rate of convergence, scored with the divergence between the current and
% final density, increases with the relative amount of solenoidal mixing.
% This is accompanied by an increase in the information length from any
% initial density to the density in the long-term future.
%
% The second section rehearses the same mechanics but in the context of
% self-organisation. To talk about self organisation it is necessary to
% separate the self from nonself by constructing a random dynamical system
% with a Markov blanket. One can then associate the conditional density
% over external states, conditioned on particular states, with a Bayesian
% belief encoded by internal states. This corresponds to the variational
% density that underwrites the free energy principle (under some
% simplifying assumptions). The marginal density over particular (i.e.,
% internal states and their blanket) states now plays the role of a
% description of the dynamics of a particle, that shows the same
% dependencies on solenoidal flow above. Namely, increasing solenoidal flow
% decreases the path integral of free energy or the rate of convergence to
% nonequilibrium steady-state. At the same time, the information length of
% paths into the future increases. The accompanying information theoretic
% measures—of the conditional density over external states and particular
% states—can be read as an extrinsic and intrinsic information geometries,
% respectively. These conjugate geometries can, in turn, be associated with
% variational free energy and thermodynamic free energy.
%
% An increase in solenoidal flow, relative to dissipative flow, goes
% hand-in-hand with the size of a particle, where random fluctuations are
% averaged away. In other words, large particles are necessarily precise
% particles that feature solenoidal flows. These flows underwrite a rapid
% convergence to nonequilibrium steady-state from any initial conditions
% that, necessarily, entail large information lengths. In short, large,
% precise particles have an itinerant aspect to their dynamics and move
% through many discernible publicity configurations from any initial
% density. This itinerancy lends precise particles and elemental kind of
% memory, in the sense that there is a greater distance between any initial
% state and future states in the space of probability densities or beliefs.
% Please see the annotated code below for further details.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: FEP_lorenz_surprise.m 8152 2021-09-13 09:17:36Z karl $



%% density dynamics - initial and final densities
%==========================================================================
% The numerical experiments in the code below leverage the Helmholtz
% decomposition of any random dynamical system that possesses a solution to
% density dynamics (i.e., the Fokker Planck equation). Specifically, in
% what follows, we specify the dynamics in terms of three components. The
% first is the log density (S) corresponding to the final or nonequilibrium
% steady-state density, a dissipative operator (G) corresponding to the
% amplitude of random fluctuations and a conservative or solenoidal
% operator (Q) that breaks detailed balance. This parameterisation or
% decomposition of dynamics enables one to specify systems with the same
% log density (c.f., a Lyapunov function) and, implicitly, characteristic
% states. However, these states can be accessed in the setting of different
% combinations of dissipative and solenoidal flow. We will be particularly
% interested in increasing the contribution of the solenoidal operator,
% relative to the dissipative operator. To enable a quasi-analytic
% treatment, we focus on functional forms based upon polynomial expansions.
% By limiting the log density to a second order function of states, the
% solution to the density dynamics becomes Gaussian. Similarly, by limiting
% the solenoidal flow to a first-order function of states, one has a fairly
% expressive model of stochastic chaos. A particular advantage of
% parameterising dynamics in this way is that conditional dependencies can
% be specified by setting appropriate second order polynomial coefficients
% to zero. These coefficients correspond to the curvature or Hessian of the
% log density, where a zero Hessian implies conditional independence. This
% device is used to construct systems with a Markov blanket with arbitrary
% (state dependent) solenoidal flows.
%--------------------------------------------------------------------------
spm_figure('GetWin','Density dynamics');
if ~nargin
    clf
    
    FEP_information_length(1,0,1)
    
    FEP_information_length(1,1,2)
    
    return
    
end

% support of phase space (for display)
%--------------------------------------------------------------------------
K    = 3;                               % second-order low density
N    = [64 1 1 64];                     % probability bins and dimensions
n    = numel(N);                        % number of dimensions
ip   = [2,3,4];                         % particular states
d    = [16 1 1 16];                              % range of support
for i = 1:n
    x{i} = linspace(-d(i),d(i),N(i));
end

% initial (Gaussian) density
%--------------------------------------------------------------------------
m0   = zeros(1,n) + 4;                  % expectation
c0   = eye(n,n);                        % covariance
S0   = spm_ness_N2Sp(m0,c0);            % polynomial parameters
P0   = spm_ness_Sp2p(S0,x);             % corresponding density

% final (Gaussian) density
%--------------------------------------------------------------------------
mT   = zeros(1,n);                      % centre on origin
H    = [ ...                            % Hessian or precision
    4 -2, 0, 0;
   -2, 4,-2, 0;
    0,-2, 4 -2;
    0, 0,-2, 4]/4;
cT   = inv(H);
cT   = cT(1:n,1:n);
ST   = spm_ness_N2Sp(mT,cT);
PT   = spm_ness_Sp2p(ST,x);

% conditional (extrinsic) density, given expected particular states
%--------------------------------------------------------------------------
[mE,cE] = spm_ness_cond(n,K,ST,ip,mT(ip));

% particular (intrinsic) density
%--------------------------------------------------------------------------
mI      = mT(ip);
cI      = cT(ip,ip);

subplot(3,2,1)
imagesc(x{4},x{1},1 - squeeze(sum(P0,[2,3])))
axis square xy
xlabel('2nd state'), ylabel('1st state')
title('Initial density','Fontsize',14)

subplot(3,2,2)
imagesc(x{4},x{1},1 - squeeze(sum(PT,[2,3])))
axis square xy
xlabel('2nd state'), ylabel('1st state')
title('NESS density','Fontsize',14)

% solenoidal coefficients
%--------------------------------------------------------------------------
[b,D,H,o] = spm_polymtx(x,K);
np = numel(ST);                            % number of polynomial coefficients
QC = [ ...                                 % constraints on solenoidal flow
    0 1 0 0;
    1 0 0 0;
    0 0 0 1;
    0 0 1 0];
QP = {};
for i = 1:n
    for j = i:n
        qp = zeros(np,1);
        if j > i && QC(i,j)
            qp(sum(o) == 0) = 1;
            qp(sum(o) == 1) = 1/64;
        else
            qp = zeros(np,1);
        end
        QP   = [QP {qp}];
    end
end
Qp    = spm_cat(QP(:));

%% Loop over different different levels of solenoidal flow
%==========================================================================
T     = 512;                                % length of timeseries
dt    = 1/64;                               % timesteps
a     = .01;                                % memory for density graphics
color = get(gca,'ColorOrder');
col   = color(ci,:);

% parameters of Helmholtz decomposition
%--------------------------------------------------------------------------
Ep.Sp = full(ST);                       % log density
Ep.Qp = full(Qp)*qi;                    % solenoidal operator
Ep.G  = eye(n,n)*gi;                    % random fluctuations

M.W   = inv(Ep.G*2);                    % precision of random fluctuations
M.K   = K;                              % order of polynomial expansion

%  Run system forwards in time to nonequilibrium steady-state density
%--------------------------------------------------------------------------
s     = m0;                             % initial state
r     = m0;                             % initial state
w     = sqrt(2*Ep.G);                   % covariance of w

Sp    = S0;                             % log density coefficients
Q     = P0;                             % initial density
FL    = [];                             % divergence or free energy
LL    = [];                             % information length
u     = [];                             % trajectory (least action)
v     = [];                             % trajectory (stochastic)
for t = 1:T
    
    % exemplar trajectories (least action)
    %======================================================================
    
    % least action
    %----------------------------------------------------------------------
    ds     = spm_NESS_gen(Ep,M,s);                % flow
    s      = s + ds*dt;                           % update
    v(:,t) = s';                                  % expected trajectory
    
    % stochastic
    %----------------------------------------------------------------------
    dr     = spm_NESS_gen(Ep,M,r);                % flow
    r      = r  + dr*dt + randn(1,n)*w*dt;        % update
    u(:,t) = r';                                 % next state
    
    
    % density dynamics
    %======================================================================
    dS = spm_NESS_ds(Sp,Ep);
    Sp = Sp + dS*dt;
    P  = spm_softmax(b*Sp);
    
    
    % divergence from steady-state (joint)
    %----------------------------------------------------------------------
    [mt,ct] = spm_ness_Sp2N(Sp,n,K);
    FL(t,1) = spm_kl_normal(mt,ct,mT,cT);
    
    % Accumulate information length (joint)
    %----------------------------------------------------------------------
    [ms,cs] = spm_ness_Sp2N(Sp - dS*dt*2,n,K);
    LL(t,1) = sqrt(spm_kl_normal(mt,ct,ms,cs));
    
    % divergence (extrinsic): conditioned on expected particular state
    %----------------------------------------------------------------------
    [me,ce] = spm_ness_cond(n,K,Sp,ip,mt(ip));
    FL(t,2) = spm_kl_normal(me,ce,mE,cE);
    
    % Accumulate information length (extrinsic)
    %----------------------------------------------------------------------
    [mr,cr] = spm_ness_cond(n,K,Sp - dS*dt*2,ip,ms(ip));
    LL(t,2) = sqrt(spm_kl_normal(me,ce,mr,cr));
    
    % divergence (intrinsic): marginal density over a particular state
    %----------------------------------------------------------------------
    mi      = mt(ip);
    ci      = ct(ip,ip);
    FL(t,3) = spm_kl_normal(mi,ci,mI,cI);
    
    % Accumulate information length (intrinsic)
    %---------------------------------------------------------------------
    mt      = mt(ip);
    ct      = ct(ip,ip);
    ms      = ms(ip);
    cs      = cs(ip,ip);
    LL(t,3) = sqrt(spm_kl_normal(mt,ct,ms,cs));
    
    
    P = reshape(P,N)*a + Q*(1 - a);
    Q = P;
    subplot(3,2,2)
    imagesc(x{4},x{1},1 - squeeze(sum(P,[2,3])))
    axis square xy
    xlabel('2nd state'), ylabel('1st state')
    title('NESS density','Fontsize',14)
    
    % overlay trajectory
    %------------------------------------------------------------------
    hold on
    plot(v(4,:),v(1,:),'Color',col)
    plot(u(4,:),u(1,:),':','Color',col), hold off
    drawnow
    
end

% plot convergence in terms of free energy and information length
%======================================================================

% information length and Kulback-Leibler divergence
%----------------------------------------------------------------------
subplot(3,2,4), hold on
t  = (1:T)*dt;
plot(t,FL(:,1),'-','Color',col)
plot(t,cumsum(LL(:,1)),'-.','Color',col)
title('Free energy and information length','Fontsize',16)
xlabel('time'),ylabel('divergence (nats)'), box off

subplot(3,2,5), hold on
t  = (1:T)*dt;
plot(t,FL(:,2),'-','Color',col)
plot(t,cumsum(LL(:,2)),'-.','Color',col)
title('Extrinsic','Fontsize',16)
xlabel('time'),ylabel('divergence (nats)'), box off

subplot(3,2,6), hold on
t  = (1:T)*dt;
plot(t,FL(:,3),'-','Color',col)
plot(t,cumsum(LL(:,3)),'-.','Color',col)
title('Intrinsic','Fontsize',16)
xlabel('time'),ylabel('divergence (nats)'), box off


% plot exemplar trajectories in state space
%----------------------------------------------------------------------
subplot(3,2,1), hold on
plot(v(4,:),v(1,:),'Color',col),
plot(u(4,:),u(1,:),':','Color',col)
title('trajectory','Fontsize',16)

subplot(3,2,3)
plot3(v(1,:),v(2,:),v(3,:),    'Color',col), hold on
plot3(u(1,:),u(2,:),u(3,:),':','Color',col)
plot3(mT(1),mT(2),mT(3),   '.','Color',col,'MarkerSize',32)
title('State-space','Fontsize',16)
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
axis square, box off

return


%% density dynamics as a function
%--------------------------------------------------------------------------
n     = 3;
syms  x [1 3]  'real'
syms  s [10 1] 'real'
syms  h [10 1] 'real'
syms  q [60 1] 'real'
syms  G [1 3]  'real'

% equations of motion as a Matlab function
%--------------------------------------------------------------------------
[b,D] = spm_polymtx({x1,x2,x3},3);
Sp    = h;
Ep.Sp = s;
Ep.Qp = q;
Ep.G  = diag(G);
dS    = spm_NESS_ds(Sp,Ep);

return


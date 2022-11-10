function FEP_information_length(gi,qi,ci,fi)
% Demonstration of density dynamics and information length
% FORMAT FEP_information_length(gi,qi,ci,fi)
%--------------------------------------------------------------------------
% gi  - scaling of (isotropic) random fluctuations; i.e., dissipative flow
% qi  - scaling of solenoidal flow; i.e., conservative flow
% ci  - colour index for plotting
% fi  - optional flag to print functional forms
% ti  - optional flag to reverse time halfway through the simulation
%__________________________________________________________________________
% This demonstration routine illustrates the key role of solenoidal flow
% (that breaks detailed balance) in optimisation and self-organisation. The
% first section shows that increasing solenoid flow leads to mixing that
% accelerates the convergence of a random dynamical system to its
% (nonequilibrium) steady-state or free energy minimum. Heuristically,
% solenoidal flow—on the level set of on objective function (here the log
% density of the said steady state)—can be regarded as searching for
% ‘points of entry’ in state space with 'steep' gradients. The key
% observation here is that the rate of convergence, scored with the
% divergence between the current and final density, increases with the
% relative amount of solenoidal mixing. This is accompanied by an increase
% in the information length from any initial density to the density in the
% long-term future.
%
% The second section rehearses the same mechanics but in the context of
% self-organisation. To talk about self organisation it is necessary to
% separate the self from nonself by constructing a random dynamical system
% with a Markov blanket. One can then associate the conditional density
% over external states, conditioned on blanket states, with a Bayesian
% belief encoded by internal states. This corresponds to the variational
% density that underwrites the free energy principle (under some
% simplifying assumptions). The marginal density over particular (i.e.,
% internal states and their blanket) states now plays the role of a
% description of the dynamics of a particle, that shows the same
% dependencies on solenoidal flow above. Namely, increasing solenoidal flow
% decreases the path integral of divergence or the rate of convergence to
% nonequilibrium steady-state. At the same time, the information length of
% paths into the future increases. The accompanying information theoretic
% measures—of the conditional density over external states and particular
% states—can be read as unfolding in extrinsic and intrinsic information
% geometries, respectively. These conjugate geometries can, in turn, be
% associated with variational free energy and thermodynamic free energy.
%
% An increase in solenoidal flow, relative to dissipative flow, goes
% hand-in-hand with the size of a particle, where random fluctuations are
% averaged away. In other words, large particles are necessarily precise
% particles that feature solenoidal flows. These flows underwrite a rapid
% convergence to nonequilibrium steady-state from any initial conditions
% that, necessarily, entail large information lengths. In short, large,
% precise particles have an itinerant aspect to their dynamics and move
% through many discernible probabilistic configurations from any initial
% density. This itinerancy lends precise particles an elemental kind of
% memory, in the sense that running the system backwards in time evinces a
% greater number of discernible belief states. This number corresponds to
% the information length, while the rate at which discernible belief states
% emerge corresponds to the information rate.
%
% Note that we can run the system forwards and backwards in time with
% impunity because the density dynamics are deterministic (as opposed to
% any stochastic path). This behaviour can be demonstrated by calling the
% current routine with an additional argument that reverses the direction
% of time halfway through a simulation.
%
% Please see the annotated code below for further details.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging



%% density dynamics - initial and final densities
%==========================================================================
% The numerical experiments in the code below leverage the Helmholtz
% decomposition of any random dynamical system that possesses a
% steady-state solution to its density dynamics (i.e., the Fokker Planck
% equation). Specifically, in what follows, we specify the dynamics in
% terms of three components. The first is the log density (S) corresponding
% to the final or nonequilibrium steady-state density, a dissipative
% operator (G) corresponding to the amplitude of random fluctuations and a
% conservative or solenoidal operator (Q) that breaks detailed balance.
% This parameterisation or decomposition of dynamics enables one to specify
% systems with the same log density (c.f., a Lyapunov function) and,
% implicitly, characteristic states. However, these states can be accessed
% in the setting of different combinations of dissipative and solenoidal
% flows. We will be particularly interested in increasing the contribution
% of solenoidal flow, relative to the dissipation due to random
% fluctuations.
%
% To enable a quasi-analytic treatment, we focus on functional forms based
% upon polynomial expansions. By limiting the log density to a second order
% function of states, the solution to the density dynamics becomes
% Gaussian. Similarly, by limiting the solenoidal flow to a first-order
% function of states, one has a fairly expressive model of nonlinear
% stochastic dynamics. A particular advantage of parameterising dynamics in
% this way is that conditional dependencies can be specified by setting
% appropriate second order polynomial coefficients to zero. These
% coefficients correspond to the curvature or Hessian of the log density,
% where a zero element in the Hessian matrix implies conditional
% independence. This device is used to construct systems with a Markov
% blanket with arbitrary (state dependent) solenoidal flows. In what
% follows, we deal with a minimal system with four states in the following
% order: external, sensory, active and internal states. The Markov blanket
% comprises the sensory and active states.
%--------------------------------------------------------------------------
if ~nargin
    
    spm_figure('GetWin','Density dynamics'); clf;
    
    % The first simulation considers a vanilla gradient descent cast in
    % terms of convergence to a target (final) density that is driven
    % purely by random fluctuations.
    %----------------------------------------------------------------------
    FEP_information_length(2,0,1)
    
    
    % By simply adding solenoidal flow, without changing the final density
    % or random fluctuations, the speed of convergence increases
    % dramatically, along with the information length. The total
    % information length shown in the middle panel (broken line) pertains
    % to the joint density over all states, while the lower panels
    % characterise the conditional density over external states, given
    % particular states and the density over particular states. By
    % construction, the internal state is conditionally independent of the
    % external state. This means that the extrinsic divergence and
    % information length corresponds to the conditional density given
    % blanket states. This is the posterior density approximated by the
    % variational density that is parameterised by internal states. In
    % virtue of these conditional independencies, one can interpret the
    % extrinsic information length as scoring the distance travelled in a
    % (Bayesian) belief space.
    %----------------------------------------------------------------------
    FEP_information_length(2,2,2)
    
    
    % Finally, we illustrate the kind of dynamics associated with
    % macroscopic (precise) particles by suppressing the amplitude of
    % random fluctuations. This results in a generalised synchrony that may
    % or may not be chaotic. In terms of the Bayesian mechanics implicit in
    % the particular partition, the distance travelled through the
    % (extrinsic) belief space continues to increase over time. In other
    % words, from an initial (probabilistic) state, the conditional (i.e.,
    % variational) density evolves through discernible probabilistic
    % configurations (i.e., Bayesian beliefs about external states that are
    % encoded by the expected internal state)
    %----------------------------------------------------------------------
    spm_figure('GetWin','Precise particles'); clf;
    FEP_information_length(1/32,2,1,1)
    
    %----------------------------------------------------------------------
    % In this final simulation, we have switched on the flag for the
    % functional forms of the operators and flows – that are the same for
    % all the simulations, which use the same functional form for the
    % solenoidal, dissipation and hessian matrices.
    
    return
    
else
    if nargin < 3
        ci = 1;
    end
end


% for any given scaling of dissipative and solenoidal flows:
%==========================================================================
% solve for density dynamics in a four dimensional system, where (by
% construction) the first and fourth states are conditionally independent,
% given the second and third (that play the role of sensory and active
% states, respectively). The images below focus on the joint distribution
% over the external (first) and internal (fourth) states. Notice that the
% density dynamics are evaluated in the space of the (polynomial)
% coefficients of the log density. This means the probability bins used for
% display can be specified arbitrarily.

% support of phase space (for display)
%--------------------------------------------------------------------------
K    = 3;                               % second-order (K - 1) system
N    = [64 1 1 64];                     % probability bins and dimensions
n    = numel(N);                        % number of dimensions
ip   = [2,3,4];                         % particular states
d    = [8 1 1 8];                       % range of support
for i = 1:n
    x{i} = linspace(-d(i),d(i),N(i));
end

% Initial (Gaussian) density
%--------------------------------------------------------------------------
m0   = zeros(1,n) + 4;                  % expectation
c0   = eye(n,n)/2;                      % covariance
S0   = spm_ness_N2Sp(m0,c0);            % polynomial parameters
P0   = spm_ness_Sp2p(S0,x);             % corresponding density

% Final (Gaussian) density
%--------------------------------------------------------------------------
% The final (i.e., nonequilibrium steady-state) density is specified in
% terms of its precision or Hessian. The Hessian is specified here to
% ensure conditional independences that underwrite a Markov blanket. In
% other words, the zero entries of the precision matrix (H) specify the
% conditional independencies of the steady-state (NESS) density.
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

% conditional (extrinsic) final density, given expected particular states
%--------------------------------------------------------------------------
[mE,cE] = spm_ness_cond(n,K,ST,ip,mT(ip));

% particular (intrinsic) final density
%--------------------------------------------------------------------------
mI      = mT(ip);
cI      = cT(ip,ip);

% illustrate initial and final densities
%--------------------------------------------------------------------------
subplot(3,2,1)
imagesc(x{4},x{1},1 - squeeze(sum(P0,[2,3])))
axis square xy
xlabel('4th state'), ylabel('1st state')
title('Initial density','Fontsize',14)

subplot(3,2,2)
imagesc(x{4},x{1},1 - squeeze(sum(PT,[2,3])))
axis square xy
xlabel('4th state'), ylabel('1st state')
title('NESS density','Fontsize',14)

% solenoidal coefficients
%--------------------------------------------------------------------------
% Next, the solenoidal operator is constructed to ensure sparse coupling,
% consistent with the specification of a Markov blanket in terms of the
% equations of motion. in this example, we consider solenoidal flows that
% are first-order functions of the states. This means the equations of
% motion are second order in the states, because the gradients of the log
% density are first-order in the states.
%--------------------------------------------------------------------------
[b,D,H,o] = spm_polymtx(x,K);            % get the orders of expansion (o)
np = size(o,2);                          % number of polynomial coefficients
Q0 = [ ...                               % 0th-order constraints
    0 1 0 0;
    1 0 0 0;
    0 0 0 1;
    0 0 1 0];
Q1 = [ ...                               % 1st-order constraints
    0 0 0 1;
    0 0 0 1;
    1 0 0 0;
    1 0 0 0];

QP = {};
for i = 1:n
    for j = i:n
        qp = zeros(np,1);
        if Q0(i,j)
            qp(sum(o) == 0) = 1;         % state-independent solenoidal flow
            qp(sum(o) == 1) = 1/64;      % state-dependent solenoidal flow      
            k = find(Q1(i,:));           % apply first-order constraints
            qp(any(o(k,:),1)) = 0;
        else
            qp = zeros(np,1);
        end
        QP   = [QP {qp}];
    end
end
Qp    = spm_cat(QP(:));

%% integrate the system, scaling the dissipative and solenoidal flows
%==========================================================================
T     = 512;                            % length of timeseries
dt    = 1/64;                           % timesteps
a     = .01;                            % memory for density graphics
color = get(gca,'ColorOrder');
col   = color(ci,:);

% parameters of Helmholtz decomposition
%--------------------------------------------------------------------------
Ep.Sp = full(ST);                       % final log density
Ep.Qp = full(Qp)*qi;                    % solenoidal operator
Ep.G  = eye(n,n)*gi;                    % random fluctuations

M.W   = inv(2*Ep.G);                    % precision of random fluctuations
M.K   = K;                              % order of polynomial expansion

%  Run system forwards in time to nonequilibrium steady-state density
%--------------------------------------------------------------------------
s     = m0;                             % initial state
r     = m0;                             % initial state
w     = sqrt(2*Ep.G);                   % covariance of fluctuations (w)

Sp    = S0;                             % initial log density coefficients
Q     = P0;                             % initial density
FL    = [];                             % divergence (c.f., free energy)
LL    = [];                             % information rate
u     = [];                             % trajectory (least action)
v     = [];                             % trajectory (stochastic)
for t = 1:T
    
    % switch direction of time halfway through, if requested
    %----------------------------------------------------------------------
    if nargin > 4 && t == T/2, dt = -dt, end
    
    
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
    u(:,t) = r';                                  % next state
    
    
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
    LL(t,1) = sqrt(2*spm_kl_normal(mt,ct,ms,cs));
    
    % divergence (extrinsic): conditioned on expected particular state
    %----------------------------------------------------------------------
    [me,ce] = spm_ness_cond(n,K,Sp,ip,mt(ip));
    FL(t,2) = spm_kl_normal(me,ce,mE,cE);
    
    % Accumulate information length (extrinsic)
    %----------------------------------------------------------------------
    [mr,cr] = spm_ness_cond(n,K,Sp - dS*dt*2,ip,ms(ip));
    LL(t,2) = sqrt(2*spm_kl_normal(me,ce,mr,cr));
    
    % divergence (intrinsic): marginal density over a particular state
    %----------------------------------------------------------------------
    mi      = mt(ip);
    ci      = ct(ip,ip);
    FL(t,3) = spm_kl_normal(mi,ci,mI,cI);
    
    % Accumulate information length (intrinsic)
    %----------------------------------------------------------------------
    mt      = mt(ip);
    ct      = ct(ip,ip);
    ms      = ms(ip);
    cs      = cs(ip,ip);
    LL(t,3) = sqrt(2*spm_kl_normal(mt,ct,ms,cs));
    
    
    % illustrate density dynamics by showing its trace over time
    %----------------------------------------------------------------------
    P = reshape(P,N)*a + Q*(1 - a);
    Q = P;
    subplot(3,2,2)
    imagesc(x{4},x{1},1 - squeeze(sum(P,[2,3])))
    axis square xy
    xlabel('4th state'), ylabel('1st state')
    title('Density dynamics','Fontsize',14)
    

    
    % overlay trajectory
    %----------------------------------------------------------------------
    hold on
    plot(v(4,:),v(1,:),'Color',col)
    plot(u(4,:),u(1,:),':','Color',col), hold off
    drawnow
    
    MOV(t) = getframe(gca);
    
end

% plot convergence in terms of divergence and information length
%==========================================================================

% information length and Kulback-Leibler divergence
%--------------------------------------------------------------------------
subplot(3,2,4), hold on
t  = (1:T)*abs(dt);
plot(t,FL(:,1),'-','Color',col)
plot(t,cumsum(LL(:,1)),'-.','Color',col)
title('Divergence and information length','Fontsize',16)
xlabel('time'),ylabel('divergence (nats)'), box off

subplot(3,2,5), hold on
plot(t,FL(:,2),'-','Color',col)
plot(t,cumsum(LL(:,2)),'-.','Color',col)
title('Extrinsic','Fontsize',16)
xlabel('time'),ylabel('divergence (nats)'), box off

subplot(3,2,6), hold on
plot(t,FL(:,3),'-','Color',col)
plot(t,cumsum(LL(:,3)),'-.','Color',col)
title('Intrinsic','Fontsize',16)
xlabel('time'),ylabel('divergence (nats)'), box off

% plot exemplar trajectories in state space
%--------------------------------------------------------------------------
subplot(3,2,1), hold on
plot(v(4,:),v(1,:),'Color',col),
plot(u(4,:),u(1,:),':','Color',col)

subplot(3,2,3)
plot3(v(1,:),v(2,:),v(4,:),    'Color',col), hold on
plot3(u(1,:),u(2,:),u(4,:),':','Color',col)
plot3(mT(1),mT(2),mT(3),   '.','Color',col,'MarkerSize',32)
title('State-space','Fontsize',16)
xlabel('1st state'), ylabel('2nd state'), zlabel('4th state')
axis square, box off

% assign movies to each graph object
%--------------------------------------------------------------------------
subplot(3,2,2)
set(gca,'Userdata',{MOV,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

% return unless functional forms are requested with a fourth argument
%--------------------------------------------------------------------------
if nargin < 4, return, end


%% functional form of flow and Helmholtz decomposition
%==========================================================================
ns   = numel(Ep.Sp);
nq   = numel(Ep.Qp);
is   = find(abs(Ep.Sp) > exp(-16));
iq   = find(abs(Ep.Qp) > exp(-16));

syms  w 'real'
syms  x [1 n] 'real'
syms  h [numel(is) 1] 'real'
syms  s [numel(is) 1] 'real'
syms  q [numel(iq) 1] 'real'

% replace sample points with symbolic variables
%--------------------------------------------------------------------------
O.X        = x;
O.W        = diag(kron(ones(n,1),w));
O.K        = K;

sSp        = zeros(ns,1,'like',h);
sEp.Sp     = zeros(ns,1,'like',h);
sEp.Qp     = zeros(nq,1,'like',q);
sSp(is)    = -s;
sEp.Sp(is) = -h;
sEp.Qp(iq) = q;
sEp.G      = inv(O.W)/2;


% evaluate flow, flow operators and Hessians for display
%--------------------------------------------------------------------------
[F,S,Q,L,H] = spm_NESS_gen(sEp,O);

% and associated gradients (D = dS/dx) and Jacobian (J)
%--------------------------------------------------------------------------
D     = zeros(n,1,'like',h);
J     = zeros(n,n,'like',h);
for i = 1:n
    D(i,1)  = diff(S,x(i));
    J(:,i)  = diff(F(:),x(i));
end

% display in symbolic and latex format
%--------------------------------------------------------------------------
F = F'
J = J
Q = reshape(cat(1,Q{:}),n,n)
S = S
D = D
L = L'
H = H

disp('F = ')
disp(latex(F)), disp(' ')
disp('J = ')
disp(latex(J)), disp(' ')
disp('S = ')
disp(latex(S)), disp(' ')
disp('Q = ')
disp(latex(Q)), disp(' ')
disp('D = ')
disp(latex(D)), disp(' ')
disp('L = ')
disp(latex(L)), disp(' ')
disp('H = ')
disp(latex(H)), disp(' ')

% Functional form of density dynamics
%--------------------------------------------------------------------------
% the information rate can be expressed as this derivative (squared)
%--------------------------------------------------------------------------
X     = cell(1,n);
for i = 1:n
    X{i} = x(i);
end
ds    = spm_NESS_ds(sSp,sEp,X)
disp('dS/dt = ')
disp(latex(ds)), disp(' ')

return




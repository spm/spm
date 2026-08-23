function FEP_biology
% This demonstration  uses an ensemble of particles with intrinsic (Lorenz
% attractor) dynamics and (Newtonian) short-range coupling. The setup is
% used to solve for dynamics among an ensemble of particles; from which a
% Markov blanket emerges, which forms a particle at the next hierarchical
% scale. These ensemble dynamics are then used to illustrate different
% perspectives; namely, those afforded by quantum, statistical, classical
% and Bayesian mechanics. A detailed description of each of these three
% treatments precedes each of the sections in the script. these
% descriptions are in the form of a figure legend, where each section is
% summarised with a figure.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% default settings
%--------------------------------------------------------------------------
rng(0)
N    = 256;                         % number of (Lorenz) oscillators
T    = 2056;                        % number of time bins
dt   = 1/32;                        % time interval

% parameters
%--------------------------------------------------------------------------
for i = 1:N

    if i < N/2

        % external
        %------------------------------------------------------------------
        P.p(i) = 1;                % partition
        P.r(i) = 0;                % Rayleigh parameter
        P.d(i) = 1;                % diffusion or random fluctuations
        P.k(i) = 1;                % rate constant

    elseif i < (N/2 + N/8)

        % sensory
        %------------------------------------------------------------------
        P.p(i) = 2;                % partition
        P.r(i) = 4;                % Rayleigh parameter
        P.d(i) = 1/32;             % diffusion or random fluctuations
        P.k(i) = 1;                % rate constant

    elseif i < (N/2 + N/8 + N/4)

        % internal
        %------------------------------------------------------------------
        P.p(i) = 3;                % partition
        P.r(i) = 32;               % Rayleigh parameter
        P.d(i) = 1/32;             % diffusion or random fluctuations
        P.k(i) = 1;                % rate constant

    else

        % active
        %------------------------------------------------------------------
        P.p(i) = 4;                % partition
        P.r(i) = 32;               % Rayleigh parameter
        P.d(i) = 1/32;             % diffusion or random fluctuations
        P.k(i) = 1;                % rate constant

    end

end

% states
%--------------------------------------------------------------------------
x.p  = randn(2,N)*4;                % microstates (position)
x.v  = zeros(2,N);                  % microstates (velocity)
x.q  = randn(3,N)/32;               % microstates (states)
u    = zeros(1,T);                  % exogenous fluctuations

% partition indices
%--------------------------------------------------------------------------
ee   = P.p == 1;
ss   = P.p == 2;
ii   = P.p == 3;
aa   = P.p == 4;

% generate an dynamics from initial conditions
%==========================================================================
spm_figure('GetWin','Markov blanket');clf
subplot(3,2,1)

% States
%--------------------------------------------------------------------------
% Q    - history of microstates (states)
% X    - history of microstates (position)
% V    - history of microstates (velocity)
%--------------------------------------------------------------------------
[Q,X,V] = spm_biology_solve(x,u,P,T,dt);

% motility
%--------------------------------------------------------------------------
subplot(3,2,2),cla
px  = squeeze(X(1,ii,1:T));
py  = squeeze(X(2,ii,1:T));
plot(px',py',':'), hold on
plot(mean(px),mean(py),'.b','MarkerSize',8)

xlabel('Position','FontSize',12)
title('Motility','FontSize',16)
axis([-1 1 -1 1]*8)
axis square
hold off



%% Jacobian
%==========================================================================
n   = 2;                              % number of physical dimensions
t   = T;                              % final time indices
P.d = spm_zeros(P.d);                 % supress fluctuations

% get Jacobian and plot
%--------------------------------------------------------------------------
J   = spm_biology_solve(Q(:,:,t),X(:,:,t),V(:,:,t),P);
i   = 2*(n*N) + (1:(3*N));
Jq  = J(i,i);

subplot(3,2,3)
spy(Jq),  title('Jacobian: macroscopic','FontSize',14)
xlabel('electrochemical states'),ylabel('states')
axis square

subplot(3,2,4)
i  = 1:(6*3);
spy(Jq(i,i)),  title('Jacobian microscopic','FontSize',14)
xlabel('electrochemical states'),ylabel('states')
axis square


%% Jacobian and Lyapunov exponents
%==========================================================================
spm_figure('GetWin','Jacobian'); clf;

% numerical Jacobian
%--------------------------------------------------------------------------
subplot(3,2,1)
spy(J), title('Jacobian (causal coupling)','FontSize',14)
xlabel('states'),ylabel('states')
axis square

% Lyapunov exponents
%--------------------------------------------------------------------------
[e,v] = eig(J);
v     = diag(v);
[~,i] = sort(real(v),'descend');
v     = v(i);
e     = e(:,i);

% plot
%--------------------------------------------------------------------------
subplot(3,2,2), plot(v,'.','MarkerSize',16)
hold on, plot([0,0],get(gca,'YLim'),':k'), hold off
hold on, plot(get(gca,'XLim'),[0 0],':k'), hold off
title('Lyapunov exponents','FontSize',14)
xlabel('real part'),ylabel('imaginary part')
axis square

i = find(real(v) < 0);
i = i(1:64);
subplot(3,2,3), bar(-1./real(v(i)))
title('Time constants','FontSize',14)
xlabel('eigenmode'),ylabel('secs')
axis square

subplot(3,2,4), bar(abs(imag(v(i))/(2*pi)))
title('Frequency','FontSize',14)
xlabel('eigenmode'),ylabel('Hz')
axis square

subplot(3,2,5), plot(real(e(:,i)))
title('Eigenmodes','FontSize',14)
xlabel('states'),ylabel('real value')
axis square, spm_axis tight

% Lyapunov exponent and Hausdorff dimension (Kaplan-Yorke conjecture)
%--------------------------------------------------------------------------
E     = [];
for t = (T - 0):T
    J = spm_biology_solve(Q(:,:,t),X(:,:,t),V(:,:,t),P);
    E(:,end + 1) = sort(eig(J),'descend','ComparisonMethod','real');
end

% Lyapunov exponent and Hausdorff dimension (Kaplan-Yorke conjecture)
%--------------------------------------------------------------------------
LE    = mean(real(E),2);
SE    = cumsum(LE);
j     = find(SE >= 0,1,'last');
HD    = j + sum(LE(1:j))/abs(LE(j + 1));
fprintf('Hausdorff Dimension: %.2f\n',HD)


%% illustrate a quantum perspective (quantum mechanics)
%==========================================================================
% This section illustrates the quantum treatment of a single state - a
% microstate from the first internal particle of the synthetic soup. The
% aim of this example is to show how one can characterise the dynamics of
% the state in terms of the Schrödinger potential and ensuing kinetic
% energy. Furthermore, this example illustrates how one can eschew the
% solution of the Schrödinger equation using the NESS lemma. Here, we will
% consider a single (micro) state in isolation by assuming its flow is a
% (linear) mixture of the marginal or expected flow under all other states
% and some fast random fluctuations. Although we know a lot about how these
% fluctuations are generated, we will treat them as stochastic and
% sufficiently fast that the only interesting behaviour can be captured by
% the Schrödinger potential. The timeseries is shown in the upper panel in
% terms of the state (solid line - arbitrarily assigned units of metres)
% and flow (dotted line). The sample distribution of states over time was
% evaluated in terms of the NESS potential using a fourth order polynomial
% fit to the negative logarithm of the sample density over 64 bins. The
% resulting estimate is shown in the left middle panel. From these, one can
% evaluate the Schrödinger potential (left lower panel). One can then solve
% the Schrödinger equation to evaluate the wave function over position in
% state-space (middle panel) and its Fourier transform, over momentum
% (lower middle panel). The corresponding ensemble densities over position
% and momentum are shown in the right panels, superimposed upon the
% corresponding sample densities. Finally, the density over momentum
% specifies the kinetic energy. To quantify this energy (and the
% Schrödinger potential) one needs the amplitude of the random fluctuations
% - or, equivalently, the reduced mass. This can be simply computed from
% the residuals of the flow having removed its expectation or marginal
% flow. This concludes a description of how the Schrödinger equation can be
% applied to characterise nonequilibrium steady-state dynamics.
% 
% However, this is not how the results in the figure are generated: they
% are derived directly from the NESS potential without solving the
% Schrödinger equation. In other words, the ensemble density is specified
% directly by the NESS potential, which means that the wave function (and
% its Fourier transform) can be specified directly from the ensemble
% density. Here, we somewhat arbitrarily split the ensemble density into a
% symmetric Gaussian component and an asymmetric (positive) residual. We
% then simply assigned the (square root of the) two components to the real
% and imaginary parts of the wave function. This complementary derivation
% of the wave function illustrates the point made in the main text; namely,
% that one can either generate the wave function directly from ensemble
% density or one can start from the Schrödinger potential and solve the
% Schrödinger equation.
%
% The reduced mass of a hydrogen atom is 9.1056 10^-31 The reduced mass of
% a water molecule for an O–H bond is approximately 1.56 ^-26 kg
%--------------------------------------------------------------------------
spm_figure('GetWin','quantum mechanics');clf

% get positions and velocities of all states
%--------------------------------------------------------------------------
bi    = find(logical(ss | aa));      % blanket  particles
ei    = find(logical(ee));           % external particles
mi    = find(logical(ii));           % internal particles

% extract timeseries and evaluate sample NESS density (assume mm)
%--------------------------------------------------------------------------
DT    = (1:1024) + 512;
t     = DT*dt;
q     = squeeze(X(2,ei(end),DT));
q     = spm_detrend(q,4)/1000;
b     = max(abs(q))*1.2;
b     = linspace(-b,b,64 + 1);
[n,b] = histcounts(q,b(:));
b(1)  = [];
db    = (b(2) - b(1));
n     = spm_conv(n,4);
n     = n(:)/sum(n)/db;

% approximate NESS potential (with a polynomial)
%--------------------------------------------------------------------------
Lx    = -log(n + exp(-16)/db);
PP    = [b.^0 b.^1 b.^2 b.^3 b.^4 b.^5 b.^6];
B     = pinv(PP)*Lx;
Vx    = @(B,b)[b.^0 b.^1 b.^2 b.^3 b.^4 b.^5 b.^6]*B;
dV    = @(B,b)[b.^0 2*b.^1 3*b.^2 4*b.^3 5*b.^4 6*b.^5]*B(2:end);
ddV   = @(B,b)[2*b.^0 6*b.^1 12*b.^2 20*b.^3 30*b.^4]*B(3:end);
p     = exp(-Vx(B,b));
p     = p(:)/sum(p)/db;

% wave function
%--------------------------------------------------------------------------
psi   = sqrt(p);

% evaluate the flow due to the potential gradients  SI UNITS
%--------------------------------------------------------------------------
h     = 6.62607004*1e-34/(2*pi);                    % m*m*kg/s
dqdt  = gradient(q,dt);                             % m/s
dVdb  = dV(B,q)/db;                                 % /m

% use the residuals to estimate the amplitude of effective fluctuations
%--------------------------------------------------------------------------
gam   = var(dqdt - dVdb*pinv(dVdb)*dqdt)/2;         % m*m/s
m     = h/(2*gam);                                  % kg

% combine to evaluate Schrödinger potential and kinetic energy
%--------------------------------------------------------------------------
VS    = h^2/(4*m)*(dV(B,b).^2/(db^2)/2 - ddV(B,b)/(db^2));
                                                    % kg.m.m/s/s (Joule)
f     = -h/(2*m)*dV(B,b)/db;                        % m/s
KE    = (m/2)*p'*f.^2;                              % kg*m*m/s/s (Joule)
str   = sprintf('Fluctuations (Kinetic energy : %-.2e Joules; %-.2e Kg)',KE,m);

% trajectory
%--------------------------------------------------------------------------
subplot(4,1,1)
plot(t,q,t,dqdt)
xlabel('Time (secs)', 'FontSize',12)
ylabel('State and flow (m and m/s)','FontSize',12)
title(str,'FontSize',16), spm_axis tight

% position functions
%--------------------------------------------------------------------------
subplot(4,3,4)
plot(b,Vx(B,b))
xlabel('State space (m)', 'FontSize',12)
ylabel('Nats','FontSize',12)
title('NESS potential','FontSize',16)

subplot(4,3,5)
plot(b,real(psi),b,imag(psi))
xlabel('State space (m)', 'FontSize',12)
ylabel('Amplitude (a.u.)','FontSize',12)
title('Wave function','FontSize',16)

subplot(4,3,6)
plot(b,p,b,n,'-.')
xlabel('State space (m)', 'FontSize',12)
ylabel('Probability density (a.u)','FontSize',12)
title('Ensemble density','FontSize',16)

% equivalent formulations for momentum (h*k)
%--------------------------------------------------------------------------
k      = h*b/db/b(end);                             % kg*m/s
dk     = k(2) - k(1);
PSI    = fft(psi)/sqrt(h)/2;
pk     = fftshift(PSI.*conj(PSI));
pk     = pk(:)/sum(pk)/dk;
nk     = histcounts(m*dqdt,k);
nk     = nk(:)/sum(nk)/dk;

% momentum functions
%--------------------------------------------------------------------------
subplot(4,3,7)
plot(b,VS)
xlabel('State space (m)', 'FontSize',12)
ylabel('Potential (kg.m.m/s/s (Joules)','FontSize',12)
title('Schrödinger potential','FontSize',16)

subplot(4,3,8)
plot(k,real(fftshift(PSI)),k,imag(fftshift(PSI)))
xlabel('Momentum (kg.m/s)', 'FontSize',12)
ylabel('Amplitude (a.u)','FontSize',12)
title('Wave function','FontSize',16)

subplot(4,3,9)
plot(k,pk,k,[0;nk],'-.')
xlabel('Momentum (kg.m/s)', 'FontSize',12)
ylabel('Probability density (a.u)','FontSize',12)
title('Ensemble density','FontSize',16)



%% illustrate the thermodynamic perspective (stochastic mechanics)
%==========================================================================
spm_figure('GetWin','stochastic mechanics'); clf
%--------------------------------------------------------------------------
% This figure below illustrates the characterisation of the synthetic soup
% (or active matter) in terms of classical (stochastic) thermodynamics. The
% analysis is fairly simple and proceeds as follows: first, any partition
% of states (e.g., internal states) can be treated as an ensemble. In other
% words, the behaviour of any one element can be treated as if it was
% responding to the same thermodynamic potential as all remaining elements.
% One can then evaluate the thermodynamic potential that best explains the
% flow of states. Given this expected or predicted flow one can then
% evaluate the variance or amplitude of random fluctuations and, equipped
% with a mobility coefficient (for each molecular species), one can
% then evaluate the temperature at any point in time. Given the temperature
% and thermodynamic potential, one can then evaluate the free energy. As
% the states approach their random dynamic attractor (i.e., nonequilibrium
% steady-state) these density functions converge and the ensemble density
% ceases to change. At this point, it becomes the NESS density. The
% implicit changes in the ensemble density over time can be characterised
% in terms of entropy production, which can be partitioned in a number of
% ways. In the example here, we have focused on the entropy dissipated by
% probability currents that, when multiplied by temperature, corresponds to
% heat dissipation.
%    To estimate the thermodynamic and ensemble potentials, we used a
% sixth-order polynomial expansion of position (and appropriate least
% squares estimators). The thermodynamic potential is that which best
% predicts the stochastic flow of states; whereas the ensemble density best
% predicts the sample density. To obtain more efficient estimators, we also
% averaged over 256 time bins at 32 consecutive intervals during the
% evolution of the system. We started at the first time bin to illustrate
% the thermodynamic correlates of self-organisation during which the
% principal Markov blanket was formed. For interest, we repeated the
% analysis for the internal states, the blanket states and external states.
% The upper row of images shows the evolution of temperature shown using a
% (hot) colour scale with a dot at the position of the particles (in two
% dimensions). The second panel shows the corresponding evolution of
% temperature in the three ensembles as a function of time. The interesting
% thing here is that the internal (blue) and blanket (red) states start off
% at about the same temperature. However, during the course of self
% organisation, the internal states slowly increase their temperature to
% become hotter than the external states (cyan). The third panel shows the
% corresponding free energy for each of the ensembles. This is most marked
% for the external and blanket states that could be thought of as spending
% their free energy to organise the internal states. This is reflected in
% the bottom panel that shows the corresponding heat dissipation, which is
% most marked for the external states, as might be guessed from the changes
% in the thermodynamic free energy. Although heat dissipation can fall to
% low levels - as nonequilibrium steady-state is approached - the
% temperature of our synthetic virus remains relatively high (here, the
% temperature reached about 300° Kelvin, which is roughly body
% temperature). This follows from the fact that random fluctuations are
% still in play - arising from intrinsic fluctuations of the internal
% states of each internal particle (at the underlying hierarchical level).
% These fluctuations disperse states over the thermodynamic potential,
% while flow down potential energy gradients reconstitutes the ensemble
% density. Effectively, this flow (times distance) produces heat that is
% endowed by the intrinsic fluctuations. At nonequilibrium steady-state
% these two processes are in balance and heat dissipation is eliminated;
% because probability currents are zero at all points in state space. The
% lower panel shows the corresponding entropy as a function of time. The
% external states increase their entropy initially and then entropy falls
% as the system finds its random dynamical attractor. Note that the entropy
% of states (that are destined to become densely coupled internal states)
% progressively falls; thereby, violating the second law. This is what we
% would expect in this far from equilibrium scenario. This sort of analysis
% provides an intuitive characterisation of stochastic dynamics in terms of
% constructs that underpin the first and second laws of thermodynamics.
%--------------------------------------------------------------------------

% get positions and velocities of all states
%--------------------------------------------------------------------------
bi    = find(logical(ss | aa));      % blanket  particles
ei    = find(logical(ee));           % external particles
mi    = find(logical(ii));           % internal particles

% evaluate surprise (NESS potential) with distribution over states and time
%--------------------------------------------------------------------------
nt    = 32;                          % number of evaluations
wt    = 8;                           % interval between evaluations
lt    = 512;                         % trajectory length

% Ion mobility Coefficient of Air: 0.0002
%-------------------------------------------------------------------------
kB    = 1.38064852e-8;               % Boltzmann constant fJ/K = pg.m.m/s/s/K
mu    = [.4, .2, .02];               % diffusion constant s/pg

% Stochastic dynamics of partitions
%--------------------------------------------------------------------------
nb    = 128;
J0    = (linspace(-1,1,nb).^32)';
ij    = {mi, bi, ei};
col   = {'b','r','c'};
for j = 1:numel(ij)
    
    % Stochastic dynamics of trajectories
    %----------------------------------------------------------------------
    S     = zeros(1,nt);
    for i = 1:nt

        % get trajectories (and stochastic flow)
        %------------------------------------------------------------------
        t     = (i - 1)*wt + (1:lt);
        q     = squeeze(X(2,ij{j},t))*1e-3;      % position (m)
        p     = gradient(q,dt);                  % motion (c.f., momentum)
        
        % stochastic density
        %------------------------------------------------------------------
        b     = max(abs(q(:)));
        b     = linspace(-b,b,nb);      
        [n,b] = hist(q(:),b(:));
        n     = n(:)/sum(n);
        db    = b(2) - b(1);
        
        % estimate potential and amplitude of fluctuations
        %------------------------------------------------------------------
        In    = -log(n + 1/512);                 % sample surprisal
        
        Vx    = @(b)[b.^0 b.^1 b.^2 b.^3 b.^4 b.^5];
        dVdx  = @(b)[0*b b.^0 2*b.^1 3*b.^2 4*b.^3 5*b.^4];
        XJ    = Vx(b);
        XV    = -mu(j)*dVdx(q(:));
        BJ    = pinv(XJ)*In;           % polynomial coefficients: surprisal
        BV    = pinv(XV)*p(:);         % polynomial coefficients: potential
        
        
        % evaluate gamma (temperature) by predicting flow
        %------------------------------------------------------------------
        G     = var(p(:) - XV*BV)/2;   % amplitude of fluctuations m.m/s
        Ts    = G/(mu(j)*kB);          % temperature K
        kT    = kB*Ts;                 % kB*Ts 

        % evaluate free energy and entropies
        %------------------------------------------------------------------
        Jj    = Vx(b)*BJ + J0;                 % ensemble potential
        U     = Vx(b)*BV + J0*kT;              % thermodynamic potential
        Jv    = U/kT;                          % scaled potential
        Z     = sum(exp(-Jv));                 % partition function
        Jj    = -log(spm_softmax(-Jj));        % ensemble potential
        Jv    = -log(spm_softmax(-Jv));        % thermodynamic potential
        Pj    = spm_softmax(-Jj);              % ensemble density over b
        Pv    = spm_softmax(-Jv);              % potential density over b
        Fs    = -gradient(U,db);               % thermodynamic force
        Fv    = -kT*log(Z);                    % steady state free energy
        f     = mu(j)*Fs;                      % predicted flow

        % evaluate entropy and heat production
        %------------------------------------------------------------------
        dJdx  = gradient(Jj,db);               % gradient of surprisal
        dPdx  = -dJdx.*Pj;                     % gradient of ensemble density
        jp    = f.*Pj - G*dPdx;                % probability current
        
        S(i)     = Pj'*(-log(Pj))*kB;          % entropy - kB*H
        Sflow(i) = f'*dPdx;                    % entropy - flow
        Sfluc(i) = G*dPdx'*(dPdx./Pj);         % entropy - fluctuations
        Stota(i) = jp'*(jp./Pj)/G;             % entropy - total
        Sdiss(i) = jp'*Fs/Ts;                  % entropy - dissipative
        Qt(i)    = jp'*Fs;                     % heat dissipation
        TS(i)    = Ts;                         % temperature
        Fz(i)    = Pj'*U - Ts*S(i);            % thermodynamic free energy
        Fz(i)    = kT*Pj'*(Jv - Jj) + Fv;      % thermodynamic free energy
        
        % heat maps
        %------------------------------------------------------------------
        Hm(j,i) = TS(i);
        xi{j,i} = squeeze(X(1,ij{j},t(lt/8)));
        xj{j,i} = squeeze(X(2,ij{j},t(lt/8)));
        
    end
 
    % plot results
    %----------------------------------------------------------------------     
    tt  = (1:length(Fz))*wt*dt;
    
    subplot(5,1,2), hold on
    plot(tt,TS - 273.15,col{j}), ylabel('Centigrade')
    title('Temperature','FontSize',16), spm_axis tight
    
    subplot(5,1,3), hold on
    plot(tt,Fz,col{j}), ylabel('fJ')
    title('Thermodynamic Free energy','FontSize',16), spm_axis tight
    
    subplot(5,1,4), hold on
    plot(tt,Qt,col{j}), ylabel('fJ/s')
    title('Heat dissipation','FontSize',16), spm_axis tight
    
    subplot(5,1,5), hold on
    plot(tt,S,'Color',col{j})
    xlabel('time (seconds)'), ylabel('fJ/K')
    title('Entropy','FontSize',16), spm_axis tight
    
end

legend({'internal','blanket','external'})

% heat maps (temperature)
%--------------------------------------------------------------------------
rgb   = colormap('hot');
nk    = 5;
kk    = fix(linspace(1,nt,nk));
for k = 1:nk
    subplot(5,nk,k)
    i     = kk(k);
    hm    = Hm(:,i);
    hm    = fix(1 + 63*hm/max(hm));
    for j = 1:size(Hm,1)
        plot(xi{j,i},xj{j,i},'.','MarkerSize',8,'Color',rgb(hm(j),:))
        hold on, axis square
        axis([-1 1 -1 1]*8)
    end
    set(gca,'Color','k')
end

%% plot results (densities and thermodynamic potentials).
%--------------------------------------------------------------------------
% The next figure provides illustrative potentials and density functions
% from the analysis above. In brief, they reflect the characterisation of
% stochastic thermodynamics of the external states of the ensemble. These
% graphs illustrate the relationships between distributions and flows that
% underpin quantification in terms of classical (stochastic)
% thermodynamics. In brief, this involves estimating two functions of phase
% or state space. The first (shown in red) is the surprise or
% self-information that characterises the ensemble density. The second
% (shown in blue) is a homologous potential energy function, whose
% gradients predict the flow at each point in phase space. At
% nonequilibrium steady-state, these two functions converge. In other
% words, the thermodynamic potential becomes self-information (to within an
% additive constant = thermodynamic free energy). The middle panel shows
% estimates of the ensemble and thermodynamic potentials. The ensemble
% potential (i.e., surprise) was estimated using a polynomial approximation
% to the (log) sample distribution over an ensemble of states. The
% corresponding density functions are shown in the left panel, where the
% sample distribution is shown as a dotted line. The right panel shows the
% gradients of the thermodynamic potential (blue line) that predicts the
% flow sampled by the simulation (dots).
%--------------------------------------------------------------------------
spm_figure('GetWin','stochastic graphs'); clf
subplot(4,3,1)
plot(b,n,'k:',b,Pj,'k--',b,Pv,'k')
xlabel('State (m)'), ylabel('(a.u)')
title('Ensemble density','FontSize',16), spm_axis tight

subplot(4,3,2)
plot(b,kT*Jj,'k--',b,U,'k')
xlabel('State (m)'), ylabel('Joules')
title('Potentials','FontSize',16), spm_axis tight

j = 1:32:spm_length(q);
subplot(4,3,3)
plot(b,f,'k',q(j),p(j),'k.','MarkerSize',1)
xlabel('State (m)'), ylabel('m/s')
title('Flow','FontSize',16), spm_axis tight


%% illustrate the Lagrangian perspective (classical mechanics)
%==========================================================================
% (Classical Mechanics): this section illustrates a treatment under a classical (Hamiltonian or Lagrangian)
% perspective. By taking ensemble averages over various quantities, one can
% suppress the influence of intrinsic (random) fluctuations, thereby
% revealing conservative, Hamiltonian dynamics mediated by solenoidal flow.
% Here, we focus on the average position of particular states and ask
% whether the associated motion conforms to classical predictions. The
% picture here is of a ball rolling around in a potential well where,
% crucially, the potential well changes with external states. In
% particular, we considered a second order polynomial approximation to the
% NESS potential of position, conditioned on the position of external
% particles. By formulating this dependency as a linear mixture of external
% positions, the gradients that produce the average motion of the particle
% can be expressed as a polynomial expansion that is quadratic in the
% particular positions and linear in the external positions. The polynomial
% coefficients can then be estimated, using least squares, to best predict
% the average motion of the particle; thereby specifying the (conditional)
% ensemble density and Hamiltonian dynamics. Heuristically, this
% corresponds to characterising the average behaviour of the particle as
% the motion of a marble (or ball) in a quadratic well (or bowl) that moves
% with the external states. The resulting behaviour can then be
% characterised in terms of the ball's mass that corresponds to the
% precision (i.e., inverse variance) of motion.
% 
% Upper left panel: this phase portrait summarises the behaviour we are
% trying to explain by plotting the position (state) against velocity
% (motion). In the absence of external perturbations, the trajectory should
% be a perfect circle. However, the external states are moving the
% potential energy well to produce more erratic, although entirely
% conservative, behaviour. Middle panel: this illustrates the potential
% well in terms of the corresponding (conditional) ensemble density shown
% over time, as a function of (a linear mixture of) external states. The
% shaded area corresponds to regions of high probability density and the
% white line shows the trajectory of the position of the Markov blanket.
% The black line is the corresponding motion of the Markov blanket. This is
% a nontrivial solution, in the sense that the external states are not
% simply moving the Markov blanket states - they are inducing Hamiltonian
% motion by moving the potential energy well. Upper right panel: the
% resulting predictions of the blanket motion account almost exactly for
% its (generalised) motion (blue dots). The red line is the corresponding
% prediction for a single particle of the Markov blanket - and illustrates
% that the states of motion only becomes the motion of states when
% intrinsic fluctuations are suppressed. In other words, each member of the
% ensemble is moving somewhat erratically and actively; however, their
% collective motion can be expressed as a nearly deterministic and
% instantaneous function of their collective positions. The motion being
% predicted is orthogonal to the positions (of internal and external
% particles) upon which the predictions are based.
% 
% Using the estimates of the NESS potential afforded by the emergence of
% conservative dynamics, one can now quantify the ensemble density for any
% given external state, over both the motion of state (left) and the state
% of motion (right). Lower left panel: this shows the marginal distribution
% over position, averaged over the trajectory shown. The marginal density
% (solid line) is based on the polynomial coefficients that best predicted
% motion, while the dotted line corresponds to the sample density. Lower
% right panel: this is the equivalent ensemble density over average motion.
% The precision of this density determines the effective mass of the
% particle. In this example, if we assume that motion is expressed in
% millimetres per second (i.e., slow motion at a macromolecular scale).
% Then the effective mass, given Planck's constant, is a about a femtogram.
% This corresponds to a lightweight virus.
%--------------------------------------------------------------------------
spm_figure('GetWin','classical mechanics');clf

% recover the expected state and velocity of the Markov blanket
%--------------------------------------------------------------------------
i    = 1:1024;
t    = length(i);
b    = find(ii);                               % internal states
e    = find(ee);                               % external states
q    = mean(squeeze(X(1,b,i)))';               % mean position
p    = mean(squeeze(V(1,b,i)))';               % mean velocity
q1   = squeeze(X(1,b(1),i));                   % first position
p1   = squeeze(V(1,b(1),i));                   % first velocity
qx   = squeeze(X(1,e,i))';                     % external states

% detrend and smooth
%--------------------------------------------------------------------------
sig  = 16;
q    = spm_conv(spm_detrend(q),sig,0);
p    = spm_conv(p,sig,0);
q1   = spm_conv(spm_detrend(q1),sig,0);
p1   = spm_conv(p1,sig,0);
w    = spm_conv(spm_detrend(qx),sig,0);

Qp   = var(p);                 % inverse mass or dispersion of flow

% estimates (external) state-dependent NESS potential
%--------------------------------------------------------------------------
clear XB X1
pw    = @(w)[w.^0 w.^1];
dVdB  = @(q,w)  [q^0*pw(w) 2*q^1];
phiB  = @(B,q,w)[q^1*pw(w)   q^2]*B;
for i = 1:t
    XB(i,:) = -Qp*dVdB( q(i),w(i,:));
    X1(i,:) = -Qp*dVdB(q1(i),w(i,:));
end

% coefficients of polynomial expansion of gradients - and plot
%--------------------------------------------------------------------------
B     = pinv(XB)*p;

subplot(3,2,1)
plot(p,q)
xlabel('Average velocity (micron/s)', 'FontSize',12)
ylabel('Average position (microns)','FontSize',12)
title('Phase portrait','FontSize',16)

subplot(3,2,2)
p0  = [min(p),max(p)]*(1 + 1/8);
plot(p,XB*B,'b.',p0,p0,'b--',p1/32,X1*B/32,'r:')
xlabel('Hamiltonian prediction', 'FontSize',12)
ylabel('State of motion','FontSize',12)
title('Conservative dynamics','FontSize',16), spm_axis tight

% recover state dependent density
%--------------------------------------------------------------------------
clear phi
i     = 1.2;
qq    = linspace(min(q)*i,max(q)*i,64);
tt    = 1:4:t;
phi   = zeros(numel(tt),64);
for i = 1:numel(tt)
    for j = 1:64
        phi(i,j) = phiB(B,qq(j),w(tt(i),:));
    end
end
pd    = spm_softmax(-phi');

subplot(3,1,2), cla
imagesc(tt*dt,qq,(1 - pd)), hold on
plot(tt*dt,q(tt),'w',tt*dt,p(tt),'k'), hold off, axis xy
xlabel('Time (s)', 'FontSize',12)
ylabel('Position','FontSize',12)
title('Velocity and (conditional density) over position','FontSize',16)

% marginal density over time
%--------------------------------------------------------------------------
subplot(3,2,5)
[nq,qq] = histcounts(q,64);
qq(end) = [];
dq  = qq(2) - qq(1);
nq  = nq/sum(nq)/dq;
pq  = mean(pd,2)/dq;

plot(qq,pq,qq,nq,':')
xlabel('Position (micron)', 'FontSize',12)
ylabel('Probability density (a.u.)','FontSize',12)
title('Marginal density over position','FontSize',16), spm_axis tight

subplot(3,2,6)
[np,pp] = histcounts(p,64);
pp(end) = [];
dp  = pp(2) - pp(1);
np  = np/sum(np)/dp;
qp  = exp(-pp.^2/2/Qp);
qp  = qp/sum(qp)/dp;

plot(pp,qp,pp,np,':')
xlabel('Motion (micron/s)', 'FontSize',12)
ylabel('Probability density (a.u.)','FontSize',12)
title('Marginal density over motion','FontSize',16), spm_axis tight

% marginal density over time (in femtograms)
%--------------------------------------------------------------------------
h     = 6.62607004*1e-34/(2*pi);                   % m*m*kg/s
h     = h * 10^(6 + 6 + 3 + 15);                   % micron * micron * fg/s
mass  = h/Qp;
text(0,0,sprintf('mass %-.04f fg',mass))


%% FEP: measurement and inference (Bayesian Mechanics)
%==========================================================================
% This demonstration uses an ensemble of particles with intrinsic (Lorentz
% attractor) dynamics and (Newtonian) short-range coupling.  The focus of
% this routine is to unpack the Bayesian perspective. The crucial aspect of
% this implicit inference (and the basis of the free energy principle) is
% the existence of a conditional synchronisation manifold, when
% conditioning internal and external states on blanket states. This
% provides the basis for a mapping between internal and external states
% that can be interpreted in terms of a probabilistic representation or
% inference.
%
% This Bayesian perspective is illustrated in terms of a mapping between
% the canonical modes of internal and external states (as approximated with
% a polynomial expansion). The canonical modes are evaluated using an
% estimate of the conditional expectations based upon the Euclidean
% proximity of blanket states. The ensuing posterior over external states
% is than illustrated, in relation to the actual external states. We also
% simulate event related potentials by identifying several points in time
% when the Markov blankets revisit the same neighbourhood.
%__________________________________________________________________________
 
% States
%--------------------------------------------------------------------------
% Q    - history of microstates (electrochemical)
% X    - history of microstates (position)
% V    - history of microstates (velocity)

% all kinds of states
%--------------------------------------------------------------------------
for i = 1:size(X,3)
    A(:,:,i) = [Q(:,:,i); X(:,:,i); V(:,:,i)];
end

% partition indices
%--------------------------------------------------------------------------
ee   = P.p == 1; ei = find(ee);
ss   = P.p == 2; si = find(ss);
ii   = P.p == 3; mi = find(ii);
aa   = P.p == 4; ai = find(aa);

% illustrate the Bayesian perspective (predictability of external states)
%==========================================================================
spm_figure('GetWin','Bayesian perspective');clf

ne = 32;
r  = 2;
for i = 1:ne
    dm(1,i) = r*cos(2*pi*i/ne);
    dm(2,i) = r*sin(2*pi*i/ne);
end

% get motion of external and electrochemical dysnmics of internal states
%--------------------------------------------------------------------------
T     = 512;                                      % length of timeseries
t     = size(X,3) - T;
rV    = zeros(size(V));
rX    = zeros(size(X));
for i = 1:T

    % spatial mixtures of external states
    %----------------------------------------------------------------------
    m(1,1) = mean(X(1,mi,i + t));                  % mean position
    m(2,1) = mean(X(2,mi,i + t));                  % mean position

    for j = 1:ne
        for k = 1:numel(ei)
            rX(:,j,i + t) = m + dm(:,j);
            d             = X(:,ei(k),i + t) - rX(:,j,i + t);
            pk(k)         = exp(-(d'*d)/2);
            rV(:,j,i + t) = rV(:,j,i) + pk(k)*V(:,k,i + t);
        end

        % normalise
        %------------------------------------------------------------------
        rV(:,j,i + t) = rV(:,j,i + t)/sum(pk);
    end
end

% extract states
%--------------------------------------------------------------------------
mi    = mi(end);
ei    = 1:ne;
V     = rV;
for i = 1:T
    Xe(i,:) = spm_vec(V(:,ei,i + t));              % external states (V)
    Xb(i,:) = spm_vec(A(:,[si,ai],i + t));         % blanket states  (A)
    Xm(i,:) = spm_vec(Q(:,mi,i + t));              % internal states (Q)
end

% probabilistic proximity in blanket space
%--------------------------------------------------------------------------
iC    = spm_inv(cov(Xb));
for i = 1:T
    for j = 1:T
        r      = Xb(i,:) - Xb(j,:);
        w(i,j) = exp(-(r*iC*r')/2);
    end
end

% convert into proper probability distribution at each time point
%--------------------------------------------------------------------------
w     = diag(sum(w,2))\w;

% and evaluate conditional expecatation of external and internal states
%--------------------------------------------------------------------------
xe    = zeros(size(Xe));
xm    = zeros(size(Xm));
for i = 1:T
    for j = 1:T
        xe(i,:) = xe(i,:) + w(i,j)*Xe(j,:);
        xm(i,:) = xm(i,:) + w(i,j)*Xm(j,:);
    end
end

% normalise and identify canonical eigenvariates
%--------------------------------------------------------------------------
xe    = spm_detrend(xe);
xm    = spm_detrend(xm);
CVA   = spm_cva(xe,xm);

% show results - canonical vectors over elements (mode M)
%--------------------------------------------------------------------------
subplot(3,2,1), cla

M     = 1;
Ve    = CVA.V(:,M);
Ve    = spm_unvec(Ve,V(:,ei,1));
ve    = sum(Ve.^2);
ve    = ve/max(ve);
for k = 1:length(Ve)
    c = [0 1 1]*ve(k) + [1 1 1]*(1 - ve(k));
    plot(rX(1,ei(k),end),rX(2,ei(k),end),'.','MarkerSize',32,'Color',c), hold on
end

% overplot mode of motion
%--------------------------------------------------------------------------
quiver(rX(1,ei,end),rX(2,ei,end),Ve(1,:),Ve(2,:))

Vm    = CVA.W(:,M);
Vm    = spm_unvec(Vm,Q(:,mi,1));
vm    = sum(Vm.^2);
vm    = vm/max(vm);
for k = 1:numel(mi)
    c = [0 0 1]*vm(k) + [1 1 1]*(1 - vm(k));
    plot(X(1,mi(k),end),X(2,mi(k),end),'.','MarkerSize',32,'Color',c), hold on
end
xlabel('Position', 'FontSize',12)
ylabel('Position','FontSize',12)
title('Canonical modes','FontSize',16)

n     = 3;
for i = 1:n
    subplot(6,n,n*3 + i), cla
    Ve    = CVA.V(:,i);
    Ve    = spm_unvec(Ve,V(:,ei,1));
    ve    = sum(Ve.^2);
    ve    = ve/max(ve);
    for k = 1:length(Ve)
        c = [0 1 1]*ve(k) + [1 1 1]*(1 - ve(k));
        plot(rX(1,ei(k),end),rX(2,ei(k),end),'.','MarkerSize',32,'Color',c), hold on
    end
    xlabel('Position', 'FontSize',12)
    ylabel('Position','FontSize',12)
end

% conditional synchronisation manifold (polynomial approximation)
%==========================================================================

% polynomial approximation
%--------------------------------------------------------------------------
xX    = xm*CVA.W(:,M);
xY    = xe*CVA.V(:,M);
XX    = [xX.^0 xX.^1 xX.^2 xX.^3];
bE    = pinv(XX)*xY;
qE    = XX*bE;

% conditional expectation and variance
%--------------------------------------------------------------------------
qC    = ones(size(qE));
qC    = abs(var(xY - qE)*qC/mean(qC));

% show results - conditional synchronisation manifold
%--------------------------------------------------------------------------
subplot(3,2,2), cla
plot(xX,xY,'.c' ), hold on
plot(xX,qE,'.b'), hold off
xlabel('Internal mode', 'FontSize',12)
ylabel('External mode','FontSize',12)
title('Synchronisation manifold','FontSize',16), spm_axis tight

% show results - conditional distributions as a function of time
%--------------------------------------------------------------------------
subplot(6,1,3), cla
spm_plot_ci(qE',qC(:)'),  hold on
plot(xY,'c' ), hold off
xlabel('Time', 'FontSize',12)
ylabel('External states','FontSize',12)
title('Inferred and real motion','FontSize',16), spm_axis tight

% numerical check on conditional independence
%--------------------------------------------------------------------------
subplot(3,2,6), cla
B     = [Xb]; %%% Xb.^2];
R     = eye(T,T) - B*pinv(B);
plot(R*CVA.w(:,1),R*Xe*CVA.V(:,1),'.b' )
xlabel('Internal mode', 'FontSize',12)
ylabel('External mode','FontSize',12)
title('Conditional independence','FontSize',16), spm_axis tight

% event related potentials
%==========================================================================

%  identify maxima using the external canonical variate
%--------------------------------------------------------------------------
u     = CVA.v(:,M);
ue    = [];
for i = 1:8
    [d,j]   = max(u(33:end - 65));
    j       = j + 32;
    k       = (j - 32):(j + 64);        % peristimulus time (around j)
    u((j - 8):(j + 8)) = -Inf;          % eliminate from next max(u(:,1))
    ue(:,i) = Xe(k,:)*CVA.V(:,M);
    um(:,i) = Xm(k,:)*CVA.W(:,M);
    us(i)   = j;
    
end
j    = any(ue);
us   = us(j);
um   = spm_detrend(um(:,j));

% plot onsets on conditional density
%--------------------------------------------------------------------------
subplot(6,1,3)
uy    = get(gca,'Ylim'); hold on
for i = 1:length(us),plot([1 1]*us(i),uy,':k'), end, hold off

% show time-locked (internal and external) fluctuations and their mean
%--------------------------------------------------------------------------
subplot(3,2,5), cla
pst   = 1000*(-32:64)*dt;
plot(pst,um,':b'),     hold on
plot(pst,mean(um,2),'b'), hold on
plot([0 0],get(gca,'YLim'),'--b'), hold off, axis square
xlabel('Time (msec)', 'FontSize',12)
ylabel('Electrochemical response','FontSize',12)
title('Simulated ERP','FontSize',16)

% canonical correlations
%--------------------------------------------------------------------------
% subplot(3,2,6)
% bar(CVA.r,1/2)
% xlabel('Mode', 'FontSize',12)
% ylabel('Correlation','FontSize',12)
% title('Canonical correlations','FontSize',16), axis square, spm_axis tight


return


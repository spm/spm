function FEP_physics
% This demonstration  uses an ensemble of particles with intrinsic (Lorenz
% attractor) dynamics and (Newtonian) short-range coupling. the setup is
% used to solve for dynamics among an ensemble of particles; from which a
% Markov blanket emerges (which forms a particle at the next hierarchical
% scale. These ensemble dynamics are then used to illustrate different
% perspectives; namely, those afforded by quantum, statistical and
% classical mechanics. A detailed description of each of these three
% treatments precedes each of the sections in the script. these
% descriptions are in the form of a figure legend, where each section is
% summarised with a figure.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: FEP_physics.m 7224 2017-11-18 18:10:09Z karl $


% default settings (GRAPHICS sets movies)
%--------------------------------------------------------------------------
rng('default')

% Demo of synchronization manifold using coupled Lorenz attractors
%==========================================================================
N    = 128;                         % number of (Lorenz) oscillators
T    = 2048;                        % number of time bins
dt   = 1/32;                        % time interval

% parameters
%--------------------------------------------------------------------------
P.k  = 1 - exp(-rand(1,N)*4);       % variations in temporal scale
P.a  = rand(1,N) > 2/3;             % no out influences
P.e  = exp(0);                      % energy parameter (well depth)
P.d  = 1/8;                         % amplitude of random fluctuations

% states
%--------------------------------------------------------------------------
x.p  = randn(2,N)*4;                % microstates (position)
x.v  = zeros(2,N);                  % microstates (velocity)
x.q  = randn(3,N)/32;               % microstates (states)
u    = zeros(1,T);                  % exogenous fluctuations


% generate an dynamics from initial conditions
%==========================================================================
spm_figure('GetWin','Markov blanket');clf
subplot(2,2,1)

% States
%--------------------------------------------------------------------------
% Q    - history of microstates (states)
% X    - history of microstates (position)
% V    - history of microstates (velocity)
%--------------------------------------------------------------------------
% [Q,X,V,A] = spm_soup(x,u,P,T,dt,1);
[Q,X,V,A] = spm_Manifold_solve(x,u,P,T,dt,1);


% Markov blanket - parents, children, and parents of children
%==========================================================================

% Adjacency matrix
%--------------------------------------------------------------------------
t     = (T - 256):T;                              % final time indices
L     = sparse(double(any(A(:,:,t),3)))';

% internal states (defined by principle eigenvector of Markov blanket)
%--------------------------------------------------------------------------
B     = double((L + L' + L*L'));
B     = B - diag(diag(B));
v     = spm_svd(B*B',1);
[v,j] = sort(abs(v(:,1)),'descend');

% get Markov blanket and divide into sensory and active states
%--------------------------------------------------------------------------
mm    = j(1:8);                                   % internal cluster
jj    = sparse(mm,1,1,N,1);                       % internal states
bb    = B*jj & (1 - jj);                          % Markov blanket
ee    = 1 - bb - jj;                              % external states
b     = find(bb);
e     = find(ee);
m     = find(jj);
s     = b(find( any(L(b,e),2)));
a     = b(find(~any(L(b,e),2)));

% adjacency matrix - with partition underneath (LL)
%--------------------------------------------------------------------------
k       = [e; s; a; m];
LL      = L;
LL(e,e) = LL(e,e) + 1/8;
LL(s,s) = LL(s,s) + 1/8;
LL(a,a) = LL(a,a) + 1/8;
LL(m,m) = LL(m,m) + 1/8;
LL      = LL(k,k);

% plot dynamics for the initial and subsequent time periods
%--------------------------------------------------------------------------
subplot(4,1,3)
r    = 1:512;
plot(r,squeeze(Q(1,e,r)),':c'), hold on
plot(r,squeeze(Q(1,m,r)),' b'), hold off
axis([r(1) r(end) -32 32])
xlabel('Time','FontSize',12)
title('Electrochemical dynamics','FontSize',16)

subplot(4,1,4)
r    = 1:T;
plot(r,squeeze(V(1,e,r)),':c'), hold on
plot(r,squeeze(V(1,m,r)),' b'), hold off
axis([r(1) r(end) -32 32])
xlabel('Time','FontSize',12)
title('Newtonian dymanics','FontSize',16)


% Markov blanket - self-assembly
%==========================================================================
subplot(2,2,2)
imagesc(1 - LL)
axis square
xlabel('Element','FontSize',12)
xlabel('Element','FontSize',12)
title('Adjacency matrix','FontSize',16)


% follow self-assembly
%--------------------------------------------------------------------------
GRAPHICS = 0;
clear M
for i = (T - 128):T

    % plot positions
    %----------------------------------------------------------------------
    subplot(2,2,2),set(gca,'color','w')

    px = ones(3,1)*X(1,:,i) + Q([1 2 3],:,i)/16;
    py = ones(3,1)*X(2,:,i) + Q([2 3 1],:,i)/16;
    plot(px,py,'.b','MarkerSize',8), hold on
    px = X(1,e,i); py = X(2,e,i);
    plot(px,py,'.c','MarkerSize',24)
    px = X(1,m,i); py = X(2,m,i);
    plot(px,py,'.b','MarkerSize',24)
    px = X(1,s,i); py = X(2,s,i);
    plot(px,py,'.m','MarkerSize',24)
    px = X(1,a,i); py = X(2,a,i);
    plot(px,py,'.r','MarkerSize',24)

    xlabel('Position','FontSize',12)
    ylabel('Position','FontSize',12)
    title('Markov Blanket','FontSize',16)
    axis([-1 1 -1 1]*8)
    axis square, hold off, drawnow

    % save
    %----------------------------------------------------------------------
    if i > (T - 128) && GRAPHICS
        M(i - T + 128) = getframe(gca);
    end

end

% set ButtonDownFcn
%--------------------------------------------------------------------------
if GRAPHICS
    h   = findobj(gca);
    set(h(1),'Userdata',{M,16})
    set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
    xlabel('Click for Movie','Color','r')
end


%% illustrate the quantum perspective (quantum mechanics)
%==========================================================================
% This section illustrates the quantum treatment of a single state – a
% microstate from the second external particle of the synthetic soup. The
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
% terms of the state (solid line – arbitrarily assigned units of metres)
% and flow (dotted line). The sample distribution of states over time was
% evaluated in terms of the NESS potential using a six order polynomial fit
% to the negative logarithm of the sample density over 64 bins. The
% resulting estimate and its derivatives are shown in the left middle
% panel. From these, Equation (5.4) specifies the Schrödinger potential
% (left lower panel). One can then solve the Schrödinger equation to
% evaluate the wave function over position in state-space (middle panel)
% and its Fourier transform, over momentum (lower middle panel). The
% corresponding ensemble densities over position and momentum are shown in
% the right panels, superimposed upon the corresponding sample densities.
% Finally, the density over momentum specifies the kinetic energy via
% Equation (5.9). Here, the kinetic (and potential) energy was 2.29 ^
% 10-33. To quantify this energy (and the Schrödinger potential) one needs
% the amplitude of the random fluctuations – or, equivalently, the reduced
% mass. This can be simply computed from the residuals of the flow having
% removed its expectation or marginal flow. The reduced mass of this
% quantum system was 5.52 ^ 10-38. This concludes a description of how the
% Schrödinger equation can be applied to characterise nonequilibrium
% steady-state dynamics. However, this is not how the results in this
% figure were generated: they were derived directly from the NESS potential
% without solving the Schrödinger equation. In other words, the ensemble
% density is specified directly by the NESS potential, which means that the
% wave function (and its Fourier transform) can be specified directly from
% the ensemble density. Here, we somewhat arbitrarily split the ensemble
% density into a symmetric Gaussian component and an asymmetric (positive)
% residual. We then simply assigned the (square root of the) two components
% to the real and imaginary parts of the wave function. This complementary
% derivation of the wave function illustrates the point made in the main
% text; namely, that one can either generate the wave function directly
% from ensemble density or one can start from the Schrödinger potential and
% solve the Schrödinger equation.
%--------------------------------------------------------------------------
spm_figure('GetWin','quantum mechanics');clf


% extract timeseries and evaluate sample NESS density
%--------------------------------------------------------------------------
t     = (1:T)*dt;
q     = squeeze(Q(1,s(2),:));
q     = spm_detrend(q);
b     = max(abs(q));
b     = linspace(-b,b,64);
[n,b] = hist(q,b(:));
db    = b(2) - b(1);
n     = n(:)/sum(n)/db;

% approximate NESS potential (with a polynomial)
%--------------------------------------------------------------------------
Vx    = -log(n + eps);
PP    = [b.^0 b.^1 b.^2 b.^3 b.^4 b.^5 b.^6];
W     = diag(spm_softmax(-Vx));
B     = (pinv(W*PP)*W*Vx);
Vx    = @(B,b)[b.^0 b.^1 b.^2 b.^3 b.^4 b.^5 b.^6]*B;
dV    = @(B,b)[b.^0 2*b.^1 3*b.^2 4*b.^3 5*b.^4 6*b.^5]*B(2:end);
ddV   = @(B,b)[2*b.^0 6*b.^1 12*b.^2 20*b.^3 30*b.^4]*B(3:end);
p     = exp(-Vx(B,b));
p     = p(:)/sum(p)/db;

% wave function in terms of symmetric and antisymmetric parts
%--------------------------------------------------------------------------
% ps    = (p + flipud(p))/2;
% pa    = (p - flipud(p))/2;
% psi   = sqrt(ps - ps/2) + sqrt(-1)*sqrt(ps/2 + pa);
ps    = exp(-b.^2/(2*(b(1)/4)^2));
ps    = ps*p(32);
pa    = p - ps;
psi   = sqrt(ps) + sqrt(-1)*sqrt(pa).*sign(b);

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
subplot(3,1,1)
plot(t,q,t,dqdt,':')
xlabel('Time (secs)', 'FontSize',12)
ylabel('State and flow (m ans m/s)','FontSize',12)
title(str,'FontSize',16), spm_axis tight


% position functions
%--------------------------------------------------------------------------
subplot(3,3,4)
plot(b,Vx(B,b),b,dV(B,b)/db,'-.',b,ddV(B,b)/db/db,':')
xlabel('State space (m)', 'FontSize',12)
ylabel('Nats (a.u, /m and /m/m)','FontSize',12)
title('NESS potential','FontSize',16), spm_axis tight

subplot(3,3,5)
plot(b,real(psi),b,imag(psi),':')
xlabel('State space (m)', 'FontSize',12)
ylabel('Amplitude (m/m)','FontSize',12)
title('Wave function','FontSize',16), spm_axis tight

subplot(3,3,6)
plot(b,n,':',b,psi.*conj(psi))
xlabel('State space (m)', 'FontSize',12)
ylabel('Probability density (a.u)','FontSize',12)
title('Ensemble density','FontSize',16), spm_axis tight


% equivalent formulations for momentum (h*k)
%--------------------------------------------------------------------------
PSI    = fft(psi)/sqrt(h)/2;
k      = h*b/db/b(end);                             % kg*m/s
[nk,k] = hist(m*dqdt,k);
dk     = k(2) - k(1);
nk     = nk(:)/sum(nk)/dk;

% momentum functions
%--------------------------------------------------------------------------
subplot(3,3,7)
plot(b,VS)
xlabel('State space (m)', 'FontSize',12)
ylabel('Potential (kg.m.m/s/s (Joules)','FontSize',12)
title('Schrödinger potential','FontSize',16), spm_axis tight

subplot(3,3,8)
plot(k,real(fftshift(PSI)),k,imag(fftshift(PSI)),':')
xlabel('Momentum (kg.m/s)', 'FontSize',12)
ylabel('Amplitude (a.u)','FontSize',12)
title('Wave function','FontSize',16), spm_axis tight

subplot(3,3,9)
plot(k,nk,':',k,fftshift(PSI.*conj(PSI)))
xlabel('Momentum (kg.m/s)', 'FontSize',12)
ylabel('Probability density (a.u)','FontSize',12)
title('Ensemble density','FontSize',16), spm_axis tight



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
% with a mobility coefficient (here, we used a mobility of 0.2), one can
% then evaluate the temperature at any point in time. Note that this is
% possible because we are associating the distribution over the states of
% the ensemble with the equivalent statistics and would have been observed
% over time. Given the temperature and thermodynamic potential, one can
% then evaluate the free energy using the KL divergence between the
% associated thermodynamic and ensemble densities. As the states approach
% their random dynamic attractor (i.e., nonequilibrium steady-state) these
% density functions converge and the ensemble density ceases to change. At
% this point, it becomes the NESS density. The implicit changes in the
% ensemble density over time can be characterised in terms of entropy
% production, which can be partitioned in a number of ways (see main text).
% In the example here, we have focused on the entropy dissipated by
% probability currents that, when multiplied by temperature, corresponds to
% heat dissipation. In this somewhat heuristic illustration, we have (for
% simplicity) focused on the position in one dimension of the internal
% states surrounded by the principal Markov blanket. This is the small
% virus like particle in Figure 4. To estimate the thermodynamic and
% ensemble potentials, we used a fourth order polynomial expansion of
% position (and appropriate least squares estimators). The thermodynamic
% potential is that which best predicts the stochastic flow of states;
% where the ensemble density best predicts the sample density. To obtain
% more efficient estimators, we also averaged over 256 time beans at 28
% consecutive intervals during the evolution of the system. We started at
% the 32nd time then to illustrate the thermodynamic, correlates of
% self-organisation during which the principal Markov blanket was formed.
% For interest, we repeated the analysis for the internal states, the
% blanket states and external states. The upper row of images shows the
% evolution of temperature shown using a (hot) colour scale with a dot at
% the position of the particles (in two dimensions). The second panel shows
% the corresponding evolution of temperature in the three ensembles as a
% function of time. The interesting thing here is that the internal (blue)
% and blanket (red) states start off at about the same temperature.
% However, during the course of self organisation, the internal states
% slowly increase their temperature to become hotter than the external
% states (cyan). In this example, the temperature of the internal states
% ended up being about twice the temperature of the blanket states. The
% third panel shows the corresponding free energy for each of the
% ensembles. The most notable thing here is that free energy decreases with
% time as the thermodynamic and ensemble potentials approach each other.
% This is most marked for the external and blanket states that could be
% thought of as spending their free energy to organise the internal states.
% After about five seconds, there is relatively little free energy left
% within the system. This is reflected in the bottom panel that shows the
% corresponding heat dissipation, which is most marked for the external
% states, as might be guessed from the changes in the thermodynamic free
% energy. Although heat dissipation can fall to low levels – as
% nonequilibrium steady-state is approached – the temperature of our
% synthetic virus remains relatively high (here, the temperature reached
% about 300° Kelvin, which is roughly body temperature). This follows from
% the fact that random fluctuations are still in play – arising from
% intrinsic fluctuations of the internal states of each internal particle
% (at the underlying hierarchical level). These fluctuations disperse
% states over the thermodynamic potential, while flow down potential energy
% gradients reconstitutes the ensemble density. Effectively, this flow
% (times distance) produces heat that is endowed by the intrinsic
% fluctuations. At nonequilibrium steady-state these two processes are in
% balance and heat dissipation is eliminated; because probability currents
% are zero at all points in state space. The lower panel shows the
% corresponding entropy as a function of time. The external states increase
% their entropy initially and then entropy falls as the system finds its
% random dynamical attractor. Note that the entropy of states (that are
% destined to become densely coupled internal states) progressively falls;
% thereby, violating the second law. This is what we would expect in this
% far from equilibrium scenario. Clearly, this is an idealised description
% of stochastic dynamics under the ensemble assumption. In reality, the
% dynamics are much more complicated and our treatment here can be
% construed as a mean field approximation, where all the interactions
% within and between the three ensembles are summarised with a common
% thermodynamic potential and ensemble density. Despite this approximation,
% this sort of analysis provides an intuitive characterisation of
% stochastic dynamics in terms of constructs that underpin the first and
% second laws of thermodynamics.
%--------------------------------------------------------------------------

% get positions and velocities of all states
%--------------------------------------------------------------------------
bi    = find(logical(bb));           % blanket  particles
ei    = find(logical(ee));           % external particles
mi    = find(logical(jj));           % internal particles


% evaluate surprise (NESS potential) with distribution over states and time
%--------------------------------------------------------------------------
nt    = 28;                          % number of evaluations
wt    = 32;                          % interval between evaluations
lt    = 256;                         % trajectory length
n0    = 32;                          % intial evaluation

% Ion mobility Coefficient of Air: 0.0002
%-------------------------------------------------------------------------
kB    = 1.38064852e-2;               % Boltzmann constant J/K = g.nm.nm/s/s/K
mu    = .2;                          % diffusion constant s/g
mu    = mu*kB;

J0    = (linspace(-1,1,64).^32)'*8;
        
% Stochastic dynamics of partitions
%--------------------------------------------------------------------------
ii    = {mi, bi, ei};
col   = {'b','r','c'};
for j = 1:numel(ii)
    
    % Stochastic dynamics of trajectories
    %----------------------------------------------------------------------
    for i = 1:nt

        % get trajectories (and stochastic flow)
        %------------------------------------------------------------------
        t     = (i - 1)*wt + (1:lt) + n0;
        q     = squeeze(X(2,ii{j},t));
        p     = gradient(q,dt);
        
        % stochastic density
        %------------------------------------------------------------------
        b     = linspace(min(q(:)),max(q(:)),64);      
        [n,b] = hist(q(:),b(:));
        n     = n + spm_hanning(64)'*sum(n)*exp(-8);
        n     = n(:)/sum(n);
        db    = b(2) - b(1);
        
        
        % estimate potential and amplitude of fluctuations
        %------------------------------------------------------------------
        In    = -log(n);
        Wj    = diag(spm_softmax(-In));
        Wv    = 1;
        
        Vx    = @(b)[b.^0 b.^1 b.^2 b.^3 b.^4];% b.^5 b.^6 b.^7 b.^8];
        dVdx  = @(b)[0*b b.^0 2*b.^1 3*b.^2 4*b.^3];% 5*b.^4 6*b.^5 7*b.^6 8*b.^7];
        XJ    = Vx(b);
        XV    = -mu*dVdx(q(:));
        BJ    = pinv(Wj*XJ)*Wj*In;    % polynomial coefficients - surprise
        BV    = pinv(Wv*XV)*Wv*p(:);  % polynomial coefficients - potential
        
        
        % evaluate gamma (temperature) by predicting flow
        %------------------------------------------------------------------
        G     = var(p(:) - XV*BV)/2;  % amplitude of fluctuations m.m/s
        Ts    = G/mu;                 % temperature K
        
        % BV  = G*BJ/mu;              % NESS solution

        % evaluate free energy and entropies
        %------------------------------------------------------------------
        Jj    = Vx(b)*BJ    + J0;
        Jv    = Vx(b)*BV/Ts + J0;
        Jj    = -log(spm_softmax(-Jj));        % ensemble potential
        Jv    = -log(spm_softmax(-Jv));        % thermodynamic potential
        Pj    = spm_softmax(-Jj)/db;           % ensemble density over b
        Pv    = spm_softmax(-Jv)/db;           % potential density over b
        f     = -G*gradient(Jv,db);            % predicted flow
        Fs    = f/mu;                          % thermodynamic force
        
        
        % evaluate entropy and heat production
        %------------------------------------------------------------------
        dJdx  = gradient(Jj,db);               % gradient of surpise
        dPdx  = -dJdx.*Pj;                     % gradient of ensemble density
        jp    = f.*Pj - G*dPdx;                % probability current
        
        Fz(i)    = Pj'*(log(Pj) - log(Pv))*db; % thermodynamic free energy
        S(i)     = Pj'*(-log(Pj))*db;
        Sflow(i) = f'*dPdx*db;                 % entropy - flow
        Sfluc(i) = G*dPdx'*(dPdx./Pj)*db;      % entropy - fluctuations
        Stota(i) = jp'*(jp./Pj)/G*db;          % entropy - total
        Sdiss(i) = jp'*Fs/Ts*db;               % entropy - dissipative
        Qt(i)    = jp'*Fs*db;                  % heat dissipation
        TS(i)    = Ts;                         % temperature
        
        
        % heat maps
        %------------------------------------------------------------------
        Hm(j,i) = TS(i);
        xi{j,i} = squeeze(X(2,ii{j},t(end)));
        xj{j,i} = squeeze(X(1,ii{j},t(end)));
        
    end
 
    % plot results
    %----------------------------------------------------------------------     
    tt  = (1:length(Fz))*wt*dt;
    
    subplot(5,1,2), hold on
    plot(tt,TS,col{j}), ylabel('Kelvin')
    title('Temperature','FontSize',16), spm_axis tight
    
    subplot(5,1,3), hold on
    plot(tt,Fz,col{j}), ylabel('Joules')
    title('Thermodynamic Free energy','FontSize',16), spm_axis tight
    
    subplot(5,1,4), hold on
    plot(tt,Qt,col{j},tt,Sflow.*TS,':','Color',col{j}),ylabel('Joules/s')
    title('Heat dissipation','FontSize',16), spm_axis tight
    
    subplot(5,1,5), hold on
    plot(tt,S,'Color',col{j})
    xlabel('Time (seconds)'), ylabel('Joules/K')
    title('Entropy','FontSize',16), spm_axis tight
    
end

% heat maps (temperature)
%--------------------------------------------------------------------------
rgb   = colormap(hot);
kk    = fix(linspace(1,nt,4));
hm    = Hm; 
% hm(:) = hm(:) - min(hm(:)) + 1;
hm(:) = ceil(64*hm/max(hm(:)));
for k = 1:4
    subplot(5,4,k)
    i = kk(k);
    for j = 1:size(Hm,1)
        plot(xi{j,i},xj{j,i},'.','MarkerSize',8,'Color',rgb(hm(j,i),:))
        hold on, axis square
        axis([-1 1 -1 1]*8)
    end
    set(gca,'Color','k')
end

% plot results(densities and thermodynamic potentials).
%--------------------------------------------------------------------------
% The next figure provides illustrative potentials and density functions
% from the analysis above. In brief, they reflect the characterisation of
% stochastic thermodynamics of the external states of the synthetic soup.
% These graphs illustrate the relationships between distributions and flows
% that underpin quantification in terms of classical (stochastic)
% thermodynamics. In brief, this involves estimating two functions of phase
% or state space. The first (shown in red) is the surprise or
% self-information that characterises the ensemble density. The second
% (shown in blue) is a homologous potential energy function, whose
% gradients predict the flow at each point in phase space. At
% nonequilibrium steady-state, these two functions are the same. In other
% words, the thermodynamic potential becomes self-information. This means
% that the distance or, more strictly speaking, divergence from
% steady-state can be quantified in terms of the KL divergence between the
% associated probability density functions (shown on the left). When these
% densities converge, the ensemble density stops changing and becomes NESS
% density. The middle panel shows estimates of the ensemble and
% thermodynamic potentials. The ensemble potential (i.e., surprise) was
% estimated using a polynomial approximation to the (log) sample
% distribution over an ensemble of states. The corresponding density
% functions are shown in the left panel, where the sample distribution is
% shown as a dotted line. The right panel shows the gradients of the
% thermodynamic potential (blue line) that predicts the flow sampled by the
% simulation (dots).
%--------------------------------------------------------------------------
spm_figure('GetWin','stochastic graphs'); clf
subplot(4,3,1)
plot(b,n/db,'r:',b,Pj,'r',b,Pv,'b')
xlabel('State (nm)'), ylabel('(a.u)')
title('Ensemble density','FontSize',16), spm_axis tight

subplot(4,3,2)
plot(b,Jj,'r',b,Jv,'b')
xlabel('State (nm)'), ylabel('Joules')
title('Potentials','FontSize',16), spm_axis tight

j = 1:32:spm_length(q);
subplot(4,3,3)
plot(b,f,'b',q(j),p(j),'b.','MarkerSize',1)
xlabel('State (nm)'), ylabel('nm/s')
title('Flow','FontSize',16), spm_axis tight


%% illustrate the Lagrangian perspective (classical mechanics)
%==========================================================================
% (Classical Mechanics): this section illustrates a treatment of our
% primordial soup under a classical (Hamiltonian or Lagrangian)
% perspective. By taking ensemble averages over various quantities, one can
% suppress the influence of intrinsic (random) fluctuations, thereby
% revealing conservative, Hamiltonian dynamics mediated by solenoidal flow.
% Here, we focus on the average position of the particles of the principal
% Markov blanket and ask whether the associated motion conforms to
% classical predictions. The picture here is of a ball rolling around in a
% potential well where, crucially, the potential well changes with external
% states. In particular, we considered a second order polynomial
% approximation to the NESS potential of the blankets position (averaged
% over its constituent particles), conditioned on the position of all
% external particles. By formulating this dependency as a linear mixture of
% external positions, the gradients that produce the average motion of the
% Markov blanket can be expressed as a polynomial expansion that is
% quadratic in the blanket positions and linear in the external positions.
% The polynomial coefficients can then be estimated, using least squares,
% to best predict the average motion of the blanket; thereby specifying the
% (conditional) ensemble density and Hamiltonian dynamics. Heuristically,
% this corresponds to characterising the average behaviour of the Markov
% blanket as the motion of a marble (or ball) in a quadratic well (or bowl)
% that moves with the external states. The resulting behaviour can then be
% characterised in terms of the ball’s mass that corresponds to the
% precision (i.e., inverse variance) of motion. Upper left panel: this
% phase portrait summarises the behaviour we are trying to explain by
% plotting the position (state) against velocity (motion). In the absence
% of external perturbations, the trajectory should be a perfect circle.
% However, it appears that the external states are moving the potential
% energy well to produce more erratic, although entirely conservative,
% behaviour. Middle panel: this illustrates the potential well in terms of
% the corresponding (conditional) ensemble density shown over time, as a
% function of (a linear mixture of) external states. The shaded area
% corresponds to regions of high probability density and the white line
% shows the trajectory of the position of the Markov blanket. The black
% line is the corresponding motion of the Markov blanket. This is a
% nontrivial solution, in the sense that the external states are not simply
% moving the Markov blanket states – they are inducing Hamiltonian motion
% by moving the potential energy well. Upper right panel: the resulting
% predictions of the blanket motion account almost exactly for its
% (generalised) motion (blue dots). The red line is the corresponding
% prediction for a single particle of the Markov blanket – and illustrates
% that the states of motion only becomes the motion of states when
% intrinsic fluctuations are suppressed. In other words, each member of the
% blanket ensemble is moving somewhat erratically and actively; however,
% their collective motion can be expressed as a nearly deterministic and
% instantaneous function of their collective positions. This is nontrivial;
% in the sense that the motion being predicted is orthogonal to the
% positions (of blanket and external particles) upon which the predictions
% are based. Using the estimates of the NESS potential afforded by the
% emergence of conservative dynamics, one can now quantify the ensemble
% density for any given external state, over both the motion of state
% (left) and the state of motion (right). Lower left panel: this shows the
% marginal distribution over position, averaged over the trajectory shown.
% The marginal density (solid line) is based on the polynomial coefficients
% that best predicted motion, while the dotted line corresponds to the
% sample density. Lower right panel: this is the equivalent ensemble
% density over average motion. The precision of this density determines the
% effective mass of the Markov blanket. In this example, if we assume that
% motion is expressed in nanometres per millisecond (i.e., slow motion at a
% macromolecular scale). Then the effective mass, given Planck's constant,
% is 136 femtograms. This corresponds to an extremely heavy virus – or a
% rather lightweight bacterium. For example, a typical E. coli would have a
% mass of 630 fg. Had we assumed that the velocity was expressed in terms
% of metres per millisecond, the mass would have been in excess of 2 tonnes
% (assuming a classical Planck's constant of unity).
%--------------------------------------------------------------------------
spm_figure('GetWin','classical mechanics');clf

% recover the expected state and velocity of the Markov blanket
%--------------------------------------------------------------------------
i    = 512:T;
t    = length(i);
b    = find(bb);                               % blanket particles
e    = find(ee);                               % external particles
q    = mean(squeeze(X(2,b,i)))';               % mean position
p    = mean(squeeze(V(2,b,i)))';               % mean velocity
q1   = squeeze(X(2,b(1),i));                   % first position
p1   = squeeze(V(2,b(1),i));                   % first velocity
qx   =     squeeze(X(1,e,i))';                 % external states
qx   = [qx squeeze(X(2,e,i))'];


% period of oscillations
%--------------------------------------------------------------------------
[~,j]= max(abs(fft(p)));
sig  = 2*t/j/2;
sig  = 32;

% find canonical loads of external states
%--------------------------------------------------------------------------
q    = spm_conv(spm_detrend(q),sig,0);
p    = spm_conv(p,sig,0);
q1   = spm_conv(spm_detrend(q1),sig,0);
p1   = spm_conv(p1,sig,0);
w    = spm_conv(spm_detrend(qx),sig,0);

Qp   = var(p);                 % inverse mass or dispersion flow

% estimates (external) state-dependent NESS potential
%--------------------------------------------------------------------------
clear XB
pw    = @(w)[w.^0 w.^1];
dVdB  = @(q,w)  [q^0*pw(w) 2*q^1];
phiB  = @(B,q,w)[q^1*pw(w)   q^2]*B;
for i = 1:t
    XB(i,:) = -Qp * dVdB( q(i),w(i,:));
    X1(i,:) = -Qp * dVdB(q1(i),w(i,:));
end

% coefficients of polynomial expansion of gradients - and plot
%--------------------------------------------------------------------------
B     = pinv(XB)*p;

subplot(3,2,1)
plot(p,q)
xlabel('Average velocity (nm/ms)', 'FontSize',12)
ylabel('Average position (nm)','FontSize',12)
title('Phase portrait','FontSize',16)

subplot(3,2,2)
p0  = [min(p),max(p)]*(1 + 1/8);
plot(p,XB*B,'b.',p0,p0,'b--',p1/32,X1*B/32,'r:')
xlabel('Hamiltonian prediction', 'FontSize',12)
ylabel('State of motion (nm/ms)','FontSize',12)
title('Conservative dynamics','FontSize',16), spm_axis tight

% recover state dependent density
%--------------------------------------------------------------------------
qq    = linspace(min(q)*1.2,max(q)*1.2,64);
tt    = 1:4:t;
for i = 1:length(tt)
    for j = 1:64
        phi(i,j) = phiB(B,qq(j),w(tt(i),:));
    end
end

pd    = spm_softmax(-phi');

subplot(3,1,2)
imagesc(tt*dt,qq,(1 - pd)), hold on
plot(tt*dt,q(tt),'w',tt*dt,p(tt),'k'), hold off,    axis xy
xlabel('Time (ms)', 'FontSize',12)
ylabel('Position','FontSize',12)
title('Conditional density','FontSize',16)

% marginal density over time
%--------------------------------------------------------------------------
subplot(3,2,5)
nq  = hist(q,qq);
dq  = qq(2) - qq(1);
nq  = nq/sum(nq)/dq;
pq  = mean(pd,2)/dq;

plot(qq,pq,qq,nq,':')
xlabel('Position (nm)', 'FontSize',12)
ylabel('Probability density (a.u.)','FontSize',12)
title('Marginal density over state','FontSize',16), spm_axis tight

subplot(3,2,6)
[np,pp]  = hist(p,64);
dp  = pp(2) - pp(1);
np  = np/sum(np)/dp;
qp  = exp(-pp.^2/2/Qp);
qp  = qp/sum(qp)/dp;

plot(pp,qp,pp,np,':')
xlabel('Motion (nm/ms)', 'FontSize',12)
ylabel('Probability density (a.u.)','FontSize',12)
title('Marginal density over motion','FontSize',16), spm_axis tight


% marginal density over time (in femtograms)
%--------------------------------------------------------------------------
h     = 6.62607004*1e-34/(2*pi);                    % m*m*kg/s
h     = h *1e33;                                    % nm * nm * fg/ms
mass  = h/Qp;
text(0,20,sprintf('mass %-.0f fg',mass))

return


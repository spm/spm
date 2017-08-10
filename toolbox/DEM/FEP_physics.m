function FEP_physics
% This demonstration  uses an ensemble of particles with intrinsic (Lorentz
% attractor) dynamics and (Newtonian) short-range coupling.  The focus of
% this routine is to unpack the Bayesian perspective. We first simulate
% dynamics  to nonequilibrium steady-state, identify the Markov blanket and
% then examine the encoding of external states by internal states; in terms
% of their expected values.
%
% The crucial aspect of this implicit inference (and the basis of the free
% energy principle) is the existence of a conditional synchronisation
% manifold, when conditioning internal and external states on the Markov
% blanket. This provides the basis for a mapping between internal and
% external states that can be interpreted in terms of a probabilistic
% representation or inference.
%
% This Bayesian perspective is illustrated in terms of a mapping between
% the canonical modes of internal and external states (as approximated
% with a polynomial expansion). The canonical modes her are evaluated
% using an estimate of the conditional expectations  based upon the
% Euclidean proximity of Markov blanket states. The ensuing posterior over
% external states is than illustrated, in relation to the actual external
% states. We also  simulate event related potentials by identifying
% several points in time when the Markov blankets revisit the same
% neighbourhood. Finally, to illustrate the underlying dynamics, the
% Jacobians or coupling among internal and external states are
% presented; using different orders of coupling (i.e., degrees of
% separation)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: FEP_physics.m 7154 2017-08-10 20:45:05Z karl $
 
 
% default settings (GRAPHICS sets movies)
%--------------------------------------------------------------------------
rng('default')
GRAPHICS = 0;
 
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
if GRAPHICS
    subplot(2,2,1)
else
    subplot(2,1,1)
end
[Q,X,V,A,x] = spm_Manifold_solve(x,u,P,T,dt,1);

% States
%--------------------------------------------------------------------------
% Q    - history of microstates (states)
% X    - history of microstates (position)
% V    - history of microstates (velocity)

for i = 1:size(X,3)
    S(:,:,i) = [Q(:,:,i);X(:,:,i);V(:,:,i)];
end
 
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
j     = j(1:8);                                   % internal cluster
jj    = sparse(j,1,1,N,1);                        % internal states
bb    = B*jj & (1 - jj);                          % Markov blanket
ee    = 1 - bb - jj;                              % external states
b     = find(bb);
e     = find(ee);
s     = b(find( any(L(b,e),2)));
a     = b(find(~any(L(b,e),2)));
m     = j;
 
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
plot(r,squeeze(Q(1,j,r)),' b'), hold off
axis([r(1) r(end) -32 32])
xlabel('Time','FontSize',12)
title('Electrochemical dynamics','FontSize',16)
 
subplot(4,1,4)
r    = 1:T;
plot(r,squeeze(V(1,e,r)),':c'), hold on
plot(r,squeeze(V(1,j,r)),' b'), hold off
axis([r(1) r(end) -32 32])
xlabel('Time','FontSize',12)
title('Newtonian dymanics','FontSize',16)
 
 
% Markov blanket - self-assembly
%==========================================================================
subplot(2,2,1)
imagesc(1 - LL)
axis square
xlabel('Element','FontSize',12)
xlabel('Element','FontSize',12)
title('Adjacency matrix','FontSize',16)

 
% follow self-assembly
%--------------------------------------------------------------------------
clear M
for i = (T - 512):T
    
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
 
% % set ButtonDownFcn
% %--------------------------------------------------------------------------
% if GRAPHICS
%     h   = findobj(gca);
%     set(h(1),'Userdata',{M,16})
%     set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
%     xlabel('Click for Movie','Color','r')
% end
%  
% % illustrate the quantum perspective (quantum mechanics)
% %==========================================================================
% spm_figure('GetWin','quantum mechanics');clf
% 
% 
% %% extract timeseries and evaluate sample NESS density
% %--------------------------------------------------------------------------
% t     = (1:T)*dt;
% q     = squeeze(Q(1,s(2),:));
% q     = spm_detrend(q);
% b     = max(abs(q));
% b     = linspace(-b,b,64);
% [n,b] = hist(q,b(:));
% db    = b(2) - b(1);
% n     = n(:)/sum(n)/db;
% 
% % approximate NESS potential (with a polynomial)
% %--------------------------------------------------------------------------
% V     = -log(n + eps);
% PP    = [b.^0 b.^1 b.^2 b.^3 b.^4 b.^5 b.^6];
% W     = diag(V < 8);
% B     = (pinv(W*PP)*W*V);
% V     = @(B,b)[b.^0 b.^1 b.^2 b.^3 b.^4 b.^5 b.^6]*B;
% dV    = @(B,b)[b.^0 2*b.^1 3*b.^2 4*b.^3 5*b.^4 6*b.^5]*B(2:end);
% ddV   = @(B,b)[2*b.^0 6*b.^1 12*b.^2 20*b.^3 30*b.^4]*B(3:end);
% p     = exp(-V(B,b));
% p     = p(:)/sum(p)/db;
% 
% % wave function in terms of symmetric and antisymmetric parts
% %--------------------------------------------------------------------------
% % ps    = (p + flipud(p))/2;
% % pa    = (p - flipud(p))/2;
% % psi   = sqrt(ps - ps/2) + sqrt(-1)*sqrt(ps/2 + pa);
% ps    = exp(-b.^2/(2*(b(1)/4)^2));
% ps    = ps*p(32);
% pa    = p - ps;
% psi   = sqrt(ps) + sqrt(-1)*sqrt(pa).*sign(b);
% 
% % evaluate the flow due to the potential gradients  SI UNITS
% %--------------------------------------------------------------------------
% h     = 6.62607004*1e-34/(2*pi);                    % m*m*kg/s
% dqdt  = gradient(q,dt);                             % m/s
% dVdb  = dV(B,q)/db;                                 % /m
% 
% % use the residuals to estimate the amplitude of effective fluctuations
% %--------------------------------------------------------------------------
% gam   = var(dqdt - dVdb*pinv(dVdb)*dqdt)/2;         % m*m/s
% m     = h/(2*gam);                                  % kg
% 
% % combine to evaluate Schrödinger potential and kinetic energy
% %--------------------------------------------------------------------------
% VS    = h^2/(4*m)*(dV(B,b).^2/(db^2)/2 - ddV(B,b)/(db^2));
%                                                     % kg*m*m/s/s (Joule)
% f     = -h/(2*m)*dV(B,b)/db;                        % m/s   
% KE    = (m/2)*p'*f.^2;                              % kg*m*m/s/s (Joule)
% str   = sprintf('Fluctuations (Kinetic energy : %-.2e Joules; %-.2e Kg)',KE,m);
% 
% % trajectory
% %--------------------------------------------------------------------------
% subplot(3,1,1)
% plot(t,q,t,dqdt,':')
% xlabel('Time (secs)', 'FontSize',12)
% ylabel('State and flow (m ans m/s)','FontSize',12)
% title(str,'FontSize',16), spm_axis tight
% 
% 
% % position functions
% %--------------------------------------------------------------------------
% subplot(3,3,4)
% plot(b,V(B,b),b,dV(B,b)/db,'-.',b,ddV(B,b)/db/db,':')
% xlabel('State space (m)', 'FontSize',12)
% ylabel('Nats (a.u, /m and /m/m)','FontSize',12)
% title('NESS potential','FontSize',16), spm_axis tight
% 
% subplot(3,3,5)
% plot(b,real(psi),b,imag(psi),':')
% xlabel('State space (m)', 'FontSize',12)
% ylabel('Amplitude (m/m)','FontSize',12)
% title('Wave function','FontSize',16), spm_axis tight
% 
% subplot(3,3,6)
% plot(b,n,':',b,psi.*conj(psi))
% xlabel('State space (m)', 'FontSize',12)
% ylabel('Probability density (a.u)','FontSize',12)
% title('Ensemble density','FontSize',16), spm_axis tight
% 
% 
% % equivalent formulations for momentum (h*k)
% %--------------------------------------------------------------------------
% PSI    = fft(psi)/sqrt(h)/2;
% k      = h*b/db/b(end);                             % kg*m/s
% [nk,k] = hist(m*dqdt,k);
% dk     = k(2) - k(1);
% nk     = nk(:)/sum(nk)/dk;
% 
% % momentum functions
% %--------------------------------------------------------------------------
% subplot(3,3,7)
% plot(b,VS)
% xlabel('State space (m)', 'FontSize',12)
% ylabel('Potential (kg.m.m/s/s (Joules)','FontSize',12)
% title('Schrödinger potential','FontSize',16), spm_axis tight
% 
% subplot(3,3,8)
% plot(k,real(fftshift(PSI)),k,imag(fftshift(PSI)),':')
% xlabel('Momentum (kg.m/s)', 'FontSize',12)
% ylabel('Amplitude (a.u)','FontSize',12)
% title('Wave function','FontSize',16), spm_axis tight
% 
% subplot(3,3,9)
% plot(k,nk,':',k,fftshift(PSI.*conj(PSI)))
% xlabel('Momentum (kg.m/s)', 'FontSize',12)
% ylabel('Probability density (a.u)','FontSize',12)
% title('Ensemble density','FontSize',16), spm_axis tight
% 
% 
% return

%% illustrate the thermodynamic perspective (stochastic mechanics)
%==========================================================================
spm_figure('GetWin','stochastic mechanics');clf






% illustrate the Lagrangian perspective (classical mechanics)
%==========================================================================
spm_figure('GetWin','classical mechanics');clf

% recover the expected state and velocity of the Markov blanket
%--------------------------------------------------------------------------
i    = 512:T;
q    = mean(squeeze(X(2,b,i)));               % mean position
p    = mean(squeeze(V(2,b,i)));               % mean velocity
qx   =     squeeze(X(1,e,i))';                     % external states
qx   = [qx squeeze(X(2,e,i))'];
qx   = [qx squeeze(V(1,e,i))'];
qx   = [qx squeeze(V(2,e,i))'];
qx   = [qx squeeze(Q(1,e,i))'];
qx   = [qx squeeze(Q(2,e,i))'];
qx   = [qx squeeze(Q(3,e,i))'];

% find canonical loads of external states
%--------------------------------------------------------------------------
sig  = 64;
q    = spm_conv(spm_detrend(q'),sig,0);
p    = spm_conv(p',sig,0);
XE   = spm_conv(spm_detrend(qx),sig,0);
YB   = spm_detrend([q,p]);
CVA  = spm_cva(YB,XE);
w    = CVA.w;
Qp   = var(p);                       % inverse mass or dispersion flow

% estimates (external) state dependent NESS potential
%--------------------------------------------------------------------------
clear XB
pw    = @(w)[w.^0 w.^1];
dVdB  = @(q,w)  [q^0*pw(w) 2*q^1 3*q^2 4*q^3 6*q^5];
phiB  = @(B,q,w)[q^1*pw(w)   q^2   q^3   q^4   q^6]*B;
t     = length(q);
for i = 1:t
    XB(i,:) = -Qp * dVdB(q(i),w(i,:));
end

% coefficients of polynomial expansion of gradients - and plot
%--------------------------------------------------------------------------
B     = pinv(XB)*p(:);

subplot(3,2,1)
plot((1:t)*dt,q,(1:t)*dt,p,':')
xlabel('Time (seconds)', 'FontSize',12)
ylabel('Average position and velocity','FontSize',12)
title('Trajectory of averages','FontSize',16), spm_axis tight

subplot(3,2,2)
plot(p,XB*B,p,p,':')
xlabel('Hamiltonian prediction', 'FontSize',12)
ylabel('Generalised motion','FontSize',12)
title('Conservative dynamics','FontSize',16), spm_axis tight

% recover state dependent density
%--------------------------------------------------------------------------
qq    = max(abs(q))*(1 + 1/8);
qq    = linspace(-qq,qq,64);
tt    = 1:32:t;
for i = 1:length(tt)
    for j = 1:64
        phi(i,j) = phiB(B,qq(j),w(i,:));
    end
end

subplot(6,1,3)
plot((1:t)*dt,w)
ylabel('Amplitude','FontSize',12)
title('External states','FontSize',16), spm_axis tight

subplot(6,1,4)
imagesc((1:t)*dt,qq,(1 - phi'))
xlabel('Time (seconds)', 'FontSize',12)
ylabel('Position','FontSize',12)
title('Potential energy','FontSize',16)

% marginal density over time
%--------------------------------------------------------------------------
pq  = spm_softmax(-spm_vec(mean(phi)));

subplot(3,2,5)
plot(qq,pq)
xlabel('Hamiltonian prediction', 'FontSize',12)
ylabel('Generalised motion','FontSize',12)
title('Conservative dynamics','FontSize',16), spm_axis tight





        [x,y] = meshgrid(-2:.2:2, -2:.2:2);
        z = x .* exp(-x.^2 - y.^2);
        [px,py] = gradient(z,.2,.2);
        contour(z), hold on, quiver(px,py), hold off

% show results - conditional synchronisation manifold
%--------------------------------------------------------------------------
subplot(3,2,2)
plot(CVA.w(:,1),Xe*CVA.V(:,1),'.c' ), hold on
plot(CVA.w(:,1),qE,'.b' ), hold off
xlabel('Internal mode', 'FontSize',12)
ylabel('External mode','FontSize',12)
title('Synchronisation manifold','FontSize',16), spm_axis tight

% show results - conditional distributions as a function of time
%--------------------------------------------------------------------------
subplot(3,1,2)
plot(Xe*CVA.V(:,1),'c' ), hold on
spm_plot_ci(qE',qC(:)'),  hold off
xlabel('Time', 'FontSize',12)
ylabel('External states','FontSize',12)
title('Inferred and real motion','FontSize',16), spm_axis tight


return
 

function FEP_fluctuations
% This demonstration routine simulates the emergence of life - as defined
% in terms of active inference - using a synthetic primordial soup. The key
% aspect of this dynamics is that there is a separation between dynamical
% states and structural states; where the dynamical states of the
% microsystem are equipped with a Lorentz attractor and the structural
% states correspond to position and velocity. The flow of structural
% states conforms to classical Newtonian mechanics. Crucially, the physical
% motion of each microsystem is coupled to its internal dynamics and vice
% versa; where the coupling among dynamics rests upon short range
% (electrochemical) forces. This means that the dependencies among the
% dynamics of each microsystem dependent on their positions. This induces a
% dependency of the systems structural integrity on its internal dynamics -
% which leads to biological self-organisation. This biological self-
% organisation is illustrated in terms of the following:
%
% i) the existence of a Markov blanket that separates internal and external
% states, where the internal states are associated with a system that
% engages in active or embodied inference.
%
% ii) emergent inference is demonstrated by showing that the internal
% states can predict the extent states, despite their separation by the
% Markov blanket.
%
% iii) this inference (encoded by the internal dynamics) is necessary to
% maintain structural integrity, as illustrated by simulated lesion
% experiments, in which the influence of various states are quenched.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: FEP_fluctuations.m 7022 2017-02-20 10:46:37Z karl $
 
 
% default settings (GRAPHICS sets movies)
%--------------------------------------------------------------------------
rng('default')
GRAPHICS = 1;
 
% Demo of synchronization manifold using coupled Lorenz attractors
%==========================================================================
d    = 8;
N    = 128;                         % number of (Lorenz) oscillators
T    = 2048;                        % number of time bins
t    = (T - 256):T;                 % final time indices
dt   = 1/32;                        % time interval
 
% parameters
%--------------------------------------------------------------------------
P.k  = 1 - exp(-rand(1,N)*4);       % variations in temporal scale
P.a  = rand(1,N) > 2/3;             % no out influences
P.e  = exp(0);                      % energy parameter (well depth)
 
 
% states
%--------------------------------------------------------------------------
x.p  = randn(2,N)*4;                % microstates (position)
x.v  = zeros(2,N);                  % microstates (velocity)
x.q  = randn(3,N)/32;               % microstates (states)
u    = zeros(1,T);                  % exogenous fluctuations
 
 
% illustrate dynamics from initial conditions
%==========================================================================
spm_figure('GetWin','Dynamics');clf
if GRAPHICS
    subplot(2,2,1)
else
    subplot(2,1,1)
end
[Q,X,V,A,x] = spm_Manifold_solve(x,u,P,T,dt,1);


% pretty graphics
%--------------------------------------------------------------------------
if GRAPHICS
    subplot(2,2,2)
    spm_Manifold_solve(x,u,P,128,1/128,3);
end

 
% Markov blanket - parents, children, and parents of children
%==========================================================================

% Adjacency matrix
%--------------------------------------------------------------------------
L     = sparse(double(any(A(:,:,t),3)));
 
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
 
% adjacency matrix - with partition underneath (LL)
%--------------------------------------------------------------------------
k       = [e; s; a; j];
LL      = L';
LL(e,e) = LL(e,e) + 1/8;
LL(s,s) = LL(s,s) + 1/8;
LL(a,a) = LL(a,a) + 1/8;
LL(j,j) = LL(j,j) + 1/8;
LL      = LL(k,k);
 
% plot dynamics for the initial and subsequent time periods
%--------------------------------------------------------------------------
subplot(4,1,3)
r    = 1:512;
plot(r,squeeze(Q(1,e,r)),':c'), hold on
plot(r,squeeze(Q(1,j,r)),' b'), hold off
axis([r(1) r(end) -32 32])
xlabel('time','FontSize',12)
title('State dymanics','FontSize',16)
 
subplot(4,1,4)
r    = 1:2048;
plot(r,squeeze(V(1,e,r)),':c'), hold on
plot(r,squeeze(V(1,j,r)),' b'), hold off
axis([r(1) r(end) -32 32])
xlabel('time','FontSize',12)
title('State dymanics','FontSize',16)
 
 
 
% Markov blanket - self-assembly
%==========================================================================
spm_figure('GetWin','Markov blanket');
 
subplot(2,2,1)
imagesc(1 - LL)
axis square
xlabel('Element','FontSize',12)
xlabel('Element','FontSize',12)
title('Adjacency matrix','FontSize',16)
clear M
 
% follow self-assembly
%--------------------------------------------------------------------------
for i = (T-512):T
    
    % plot positions
    %----------------------------------------------------------------------
    subplot(2,2,2),set(gca,'color','w')
    
    px = ones(3,1)*X(1,:,i) + Q([1 2 3],:,i)/16;
    py = ones(3,1)*X(2,:,i) + Q([2 3 1],:,i)/16;
    plot(px,py,'.b','MarkerSize',8), hold on
    px = X(1,e,i);
    py = X(2,e,i);
    plot(px,py,'.c','MarkerSize',24)
    px = X(1,j,i);
    py = X(2,j,i);
    plot(px,py,'.b','MarkerSize',24)
    px = X(1,s,i);
    py = X(2,s,i);
    plot(px,py,'.m','MarkerSize',24)
    px = X(1,a,i);
    py = X(2,a,i);
    plot(px,py,'.r','MarkerSize',24)
    
    xlabel('Position','FontSize',12)
    ylabel('Position','FontSize',12)
    title('Markov Blanket','FontSize',16)
    axis([-1 1 -1 1]*d)
    axis square
    hold off
    drawnow
    
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
 

 
% illustrate the Bayesian perspective (predictability of external states)
%==========================================================================
spm_figure('GetWin','Bayesian perspective');clf

% establish a statistical dependency between internal (dynamic) states (XQ)
%--------------------------------------------------------------------------
m     = j;
T     = 512;                         % length of timeseries
t     = size(X,3) - T - 2;
for i = 1:T
    Xe(i,:) = spm_vec(X(:,e,i + t));
    Xb(i,:) = spm_vec(X(:,[a;s],i + t));
    Xm(i,:) = spm_vec(X(:,m,i + t));
end

xe    = zeros(size(Xe));
xE    = zeros(size(Xe));
xm    = zeros(size(Xm));
iC    = inv(cov(Xb));
for i = 1:T
    for j = 1:T
        r       = Xb(i,:) - Xb(j,:);
        w       = exp(-(r*iC*r')/16);
        xE(i,:) = xE(i,:) + w*Xe(i,:);
        xe(i,:) = xe(i,:) + w*Xe(j,:);
        xm(i,:) = xm(i,:) + w*Xm(j,:);
    end
end


% normalise and retain principal eigenvariates
%--------------------------------------------------------------------------
xE   = spm_detrend(xE);
xe   = spm_detrend(xe);
xm   = spm_detrend(xm);
CVA  = spm_cva(xe,xm);

% show results
%--------------------------------------------------------------------------
subplot(3,2,1)
V     = CVA.V(:,1);
V     = spm_unvec(V,X(:,e,1));
V     = sum(V.^2);
V     = V/max(V);
for k = 1:length(V)
    c = [0 1 1]*V(k) + [1 1 1]*(1 - V(k));
    plot(X(1,e(k),end),X(2,e(k),end),'.','MarkerSize',32,'Color',c), hold on
end

V     = CVA.W(:,1);
V     = spm_unvec(V,X(:,m,1));
V     = sum(V.^2);
V     = V/max(V);
for k = 1:length(V)
    c = [1 0 0]*V(k) + [1 1 1]*(1 - V(k));
    plot(X(1,m(k),end),X(2,m(k),end),'.','MarkerSize',32,'Color',c), hold on
end


for i = 1:T
    Xv(i) = xE(i,:)*CVA.V(:,1);
end

subplot(3,2,2)
plot(CVA.w(:,1),Xv,'.r' ),         hold on
plot(CVA.w(:,1),CVA.v(:,1),'.b' ), hold off
xlabel('Time', 'FontSize',12)
ylabel('External states','FontSize',12)
title('Motion (where)','FontSize',16)
spm_axis tight


subplot(3,1,2)
plot(CVA.v(:,1),'b:'), hold on
plot(CVA.w(:,1),'b' ), hold off
xlabel('Time', 'FontSize',12)
ylabel('External states','FontSize',12)
title('Motion (where)','FontSize',16)
spm_axis tight


return
 

function FEP_fluctuations
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
% Jacobian is or coupling among internal and external states are
% presented; using different orders of coupling (i.e., degrees of
% separation)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: FEP_fluctuations.m 7033 2017-03-05 11:19:18Z karl $
 
 
% default settings (GRAPHICS sets movies)
%--------------------------------------------------------------------------
rng('default')
GRAPHICS = 1;
 
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
P.d  = 0;                           % zero amplitude random fluctuations
 
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

 
% Markov blanket - parents, children, and parents of children
%==========================================================================

% Adjacency matrix
%--------------------------------------------------------------------------
t     = (T - 256):T;                 % final time indices
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
xlabel('time','FontSize',12)
title('State dymanics','FontSize',16)
 
subplot(4,1,4)
r    = 1:T;
plot(r,squeeze(V(1,e,r)),':c'), hold on
plot(r,squeeze(V(1,j,r)),' b'), hold off
axis([r(1) r(end) -32 32])
xlabel('time','FontSize',12)
title('State dymanics','FontSize',16)
 
 
% Markov blanket - self-assembly
%==========================================================================
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
    px = X(1,e,i); py = X(2,e,i);
    plot(px,py,'.c','MarkerSize',24)
    px = X(1,j,i); py = X(2,j,i);
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
 

% illustrate the Bayesian perspective (predictability of external states)
%==========================================================================
spm_figure('GetWin','Bayesian perspective');clf

% establish a statistical dependency between internal (dynamic) states (XQ)
%--------------------------------------------------------------------------
T     = 512;                                   % length of timeseries
t     = size(X,3) - T - 2;
for i = 1:T
    Xe(i,:) = spm_vec(X(:,e,i + t));           % external states
    Xb(i,:) = spm_vec(X(:,[a;s],i + t));       % Markov blanket
    Xm(i,:) = spm_vec(X(:,m,i + t));           % internal states
end
xe    = zeros(size(Xe));
xm    = zeros(size(Xm));
iC    = inv(cov(Xb));

% probabilistic proximity in the space of the Markov blanket
%--------------------------------------------------------------------------
for i = 1:T
    for j = 1:T
        r      = Xb(i,:) - Xb(j,:);
        w(i,j) = exp(-(r*iC*r')/32);
    end
end
w = diag(sum(w,2))\w;

% mean
%--------------------------------------------------------------------------
for i = 1:T
    for j = 1:T
        xe(i,:) = xe(i,:) + w(i,j)*Xe(j,:);
        xm(i,:) = xm(i,:) + w(i,j)*Xm(j,:);
    end
end

% covariance (not used)
%--------------------------------------------------------------------------
% ce    = zeros(size(Xe,1),size(Xe,2),size(Xe,2));
% for i = 1:T
%     for j = 1:T
%         ce(i,:,:) = squeeze(ce(i,:,:)) + w(i,j)*(Xe(i,:) - xe(i,:))'*(Xe(j,:) - xe(j,:));
%     end
% end

% normalise and identify canonical eigenvariates
%--------------------------------------------------------------------------
xe    = spm_detrend(xe);
xm    = spm_detrend(xm);
CVA   = spm_cva(xe,xm);

% show results - canonical vectors over elements
%--------------------------------------------------------------------------
subplot(3,2,1)
Ve     = CVA.V(:,1);
Ve     = spm_unvec(Ve,X(:,e,1));
Ve     = sum(Ve.^2);
Ve     = Ve/max(Ve);
for k = 1:length(Ve)
    c = [0 1 1]*Ve(k) + [1 1 1]*(1 - Ve(k));
    plot(X(1,e(k),end),X(2,e(k),end),'.','MarkerSize',32,'Color',c), hold on
end

Vm     = CVA.W(:,1);
Vm     = spm_unvec(Vm,X(:,m,1));
Vm     = sum(Vm.^2);
Vm     = Vm/max(Vm);
for k = 1:length(Vm)
    c = [1 0 0]*Vm(k) + [1 1 1]*(1 - Vm(k));
    plot(X(1,m(k),end),X(2,m(k),end),'.','MarkerSize',32,'Color',c), hold on
end


% conditional synchronisation manifold (polynomial approximation)
%==========================================================================

% polynomial approximation
%--------------------------------------------------------------------------
xX    = xm*CVA.W(:,1);
XX    = [xX.^0 xX.^1 xX.^2 xX.^3];
bE    = pinv(XX)*Xe*CVA.V(:,1);
qE    = XX*bE;

% conditional expectation and variance
%--------------------------------------------------------------------------
% for i = 1:T
%     qC(i) = CVA.V(:,1)'*squeeze(ce(i,:,:))*CVA.V(:,1);
% end
qC    = ones(size(qE));
qC    = abs(var(Xe*CVA.V(:,1) - qE)*qC/mean(qC));

% show results - conditional synchronisation manifold
%--------------------------------------------------------------------------
subplot(3,2,2)
plot(CVA.w(:,1),Xe*CVA.V(:,1),'.r' ), hold on
plot(CVA.w(:,1),qE,'.b' ), hold off
xlabel('Internal mode', 'FontSize',12)
ylabel('External mode','FontSize',12)
title('Motion (where)','FontSize',16)
spm_axis tight

% show results - conditional distributions as a function of time
%--------------------------------------------------------------------------
subplot(3,1,2)
spm_plot_ci(qE',qC(:)'),      hold on
plot(Xe*CVA.V(:,1),'r' ), hold off
xlabel('Time', 'FontSize',12)
ylabel('External states','FontSize',12)
title('Motion (where)','FontSize',16)
spm_axis tight


%  event related potentials
%==========================================================================
% for i = 1:T
%     for j = 1:T
%         r      = Xb(i,:) - Xb(j,:);
%         w(i,j) = exp(-(r*iC*r')/32);
%     end
% end

%  identify points of interest using the principal eigenvariate of w
%--------------------------------------------------------------------------
u     =  spm_svd(w);
for i = 1:8
    [d,j] = max(u(:,1));
    try
        k    = j-16:j+32;
        u(k) = -Inf;
        ue(:,i) = Xe(k,:)*CVA.V(:,1);
        um(:,i) = Xm(k,:)*CVA.W(:,1);
        us(i)   = j;
    catch
        u(j) = -Inf;
    end
end
j    = any(ue);
us   = us(j);
ue   = spm_detrend(ue(:,j));
um   = spm_detrend(um(:,j));

% plot points of interest on conditional density
%--------------------------------------------------------------------------
subplot(3,1,2)
hold on
for i = 1:length(us),plot([1 1]*us(i),[-1 1],':'),end
hold off

%  show time locked (internal and external) fluctuations and their mean
%--------------------------------------------------------------------------
subplot(3,1,3)
plot(mean(ue,2),'r' ), hold on
plot(mean(um,2),'b' ), hold on
plot(ue,'r:' ), hold on
plot(um,'b:' ), hold off
axis square


% Jacobian's and generalised synchronisation
%==========================================================================
spm_figure('GetWin','Jacobians');clf

% get Markov blanket indices the Jacobian
%--------------------------------------------------------------------------
xi  = spm_zeros(x); xi.p(:,e) = 1; iXe = find(spm_vec(xi));
xi  = spm_zeros(x); xi.p(:,m) = 1; iXm = find(spm_vec(xi));

% show results - Jacobians (of increasing order)
%--------------------------------------------------------------------------
j     = t + (1:16);
J     = spm_Manifold_solve(Q(:,:,j),X(:,:,j),V(:,:,j),P);
J     = mean(J,3);
j     = [iXe;iXm];
% q     = 8;
% U     = blkdiag(CVA.V(:,1:q),CVA.W(:,1:q));  % eigenmodes (not used)

n     = 4;
for i = 1:n
    
    % all states
    %----------------------------------------------------------------------
    JJ    = J^i;
    subplot(n,2,(i - 1)*2 + 1),spy(abs(JJ) > 1e-8)
    title(sprintf('%i-order coupling',i),'FontSize',16)
    xlabel('All states','FontSize',12)
    ylabel('All states','FontSize',12)
    
    % Internal and external states
    %----------------------------------------------------------------------
    JJ     = JJ(j,j); % JJ     = pinv(U)*JJ(j,j)*U;
    subplot(n,2,(i - 1)*2 + 2);
    spy(abs(JJ) > 1e-8)
    title(sprintf('%i-order coupling',i),'FontSize',16)
    xlabel('Ext. and int. states','FontSize',12)
    ylabel('Ext. and int.','FontSize',12)
    
end


return
 

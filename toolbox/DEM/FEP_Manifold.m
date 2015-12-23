function FEP_Manifold
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
% $Id: FEP_Manifold.m 6655 2015-12-23 20:21:27Z karl $
 
 
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
s     = b(find( any(L(e,b))));
a     = b(find(~any(L(e,b))));
 
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
 
 
% overlay attributes
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(X(1,e,T),X(2,e,T),'.c','MarkerSize',24), hold on
plot(X(1,j,T),X(2,j,T),'.b','MarkerSize',24), hold on
plot(X(1,s,T),X(2,s,T),'.m','MarkerSize',24), hold on
plot(X(1,a,T),X(2,a,T),'.r','MarkerSize',24), hold on
px = X(1,find(~P.a),T);
py = X(2,find(~P.a),T);
plot(px,py,'.w','MarkerSize',22)
xlabel('Position','FontSize',12)
title('Dynamically closed','FontSize',16)
axis([-1 1 -1 1]*d)
axis square
hold off
 
subplot(2,2,4)
plot(X(1,e,T),X(2,e,T),'.c','MarkerSize',24), hold on
plot(X(1,j,T),X(2,j,T),'.b','MarkerSize',24), hold on
plot(X(1,s,T),X(2,s,T),'.m','MarkerSize',24), hold on
plot(X(1,a,T),X(2,a,T),'.r','MarkerSize',24), hold on
px = X(1,find(P.k > 1/2),T);
py = X(2,find(P.k > 1/2),T);
plot(px,py,'.w','MarkerSize',22)
xlabel('Position','FontSize',12)
title('Slow systems','FontSize',16)
axis([-1 1 -1 1]*d)
axis square
hold off
drawnow
 
 
% Evolution of the Markov blanket
%==========================================================================
spm_figure('GetWin','Evolution of the Markov blanket');clf
rng('default')
 
T     = 512;
for i = 1:4
    
    % intervention
    %--------------------------------------------------------------------------
    PP  = P;
    if i == 2
        PP.a(a) = 1;               % lesion active elements
    elseif i == 3
        PP.a(s) = 1;               % lesion sensory elements
    elseif i == 4
        PP.a(j) = 1;               % lesion internal elements
    end
    
    % continue integration
    %--------------------------------------------------------------------------
    [QQ,XX] = spm_Manifold_solve(x,u,PP,T,dt,0);
    
    % plot trajectories of Markov Blanket and final internal state
    %--------------------------------------------------------------------------
    subplot(2,2,i),set(gca,'color','w')
    px = squeeze(XX(1,s,:));
    py = squeeze(XX(2,s,:));
    plot(px,py,'.m','MarkerSize',4), hold on
    px = squeeze(XX(1,a,:));
    py = squeeze(XX(2,a,:));
    plot(px,py,'.r','MarkerSize',4), hold on
    px = squeeze(XX(1,j,:));
    py = squeeze(XX(2,j,:));
    plot(px,py,'.b','MarkerSize',8), hold off
    xlabel('Position','FontSize',12)
    axis([-1 1 -1 1]*d)
    axis square
    
    % title
    %--------------------------------------------------------------------------
    if i == 1
        title('Markov blanket','FontSize',16)
    elseif i == 2
        title('Active lesion','FontSize',16);
    elseif i == 3
        title('Sensory lesion','FontSize',16);
    elseif i == 4
        title('Internal lesion','FontSize',16);
    end
    
end
 
 
% illustrate the Bayesian perspective (predictability of external states)
%==========================================================================
spm_figure('GetWin','Bayesian perspective');clf

% establish a statistical dependency between internal (dynamic) states (XQ)
%--------------------------------------------------------------------------
clear XQ
T     = 512;                         % length of timeseries
n     = 64;                          % temporal embedding
t     = size(Q,3) - n - T - n;
for i = 1:T
    for k = 1:(n + n)
        XQ{i,k} = spm_vec(Q(:,j,i + k + t))';
    end
end

% normalise and retain principal eigenvariates
%--------------------------------------------------------------------------
Xq        = spm_en(diff(spm_cat(XQ),1));
[Xq,S,UX] = spm_svd(Xq,1);
Xq        = Xq(:,1:32);
X0        = ones(size(Xq,1),1);

subplot(3,2,1)
imagesc(Xq')
xlabel('Time', 'FontSize',12)
ylabel('Modes','FontSize',12)
title('Internal states','FontSize',16)
axis square

% get equivalent external states to be predicted
%--------------------------------------------------------------------------
t     = size(Q,3) + (1:T) - n - T;
for c = 1:length(e)
    
    % get the external states for this element
    %----------------------------------------------------------------------
    Y    = squeeze(X(:,e(c),t))';
    Y0   = flipud(Y);
    Y    = diff(spm_en(Y ),1);
    Y0   = diff(spm_en(Y0),1);
    
    % and take the principal eigenvariates
    %----------------------------------------------------------------------
    Yq   = spm_svd(Y,1);
    Yq0  = spm_svd(Y0,1);
    
    % test for dependencies using canonical variance analysis
    %----------------------------------------------------------------------
    CVA(c)  = spm_cva(Yq,Xq,X0);
    CVA0    = spm_cva(Yq0,Xq,X0);
    SPM(c)  = CVA(c).chi(1);
    SPM0(c) = CVA0.chi(1);
    
end

% show results
%--------------------------------------------------------------------------
subplot(3,2,3)
[k,i] = max(SPM);
CVA   = CVA(i);
plot(CVA.v(:,1),'b:'), hold on
plot(CVA.w(:,1),'b' ), hold off
xlabel('Time', 'FontSize',12)
ylabel('External states','FontSize',12)
title('Motion (where)','FontSize',16)
axis square
spm_axis tight

% SPM of external prediciability
%--------------------------------------------------------------------------
subplot(3,2,4)
W     = (SPM/max(SPM)).^2;
for k = 1:length(e)
    c = [0 1 1]*W(k) + [1 1 1]*(1 - W(k));
    plot(X(1,e(k),end),X(2,e(k),end),'.','MarkerSize',32,'Color',c), hold on
end
plot(X(1,e(i),end),X(2,e(i),end),'.m','MarkerSize',32), hold on
plot(X(1,j,end),X(2,j,end),'.b','MarkerSize',32), hold off
xlabel('Position','FontSize',12)
ylabel('Position','FontSize',12)
title('Predictability','FontSize',16)
axis([-1 1 -1 1]*d)
axis square
hold off

% Distributions
%--------------------------------------------------------------------------
subplot(3,2,2)
hist([SPM' SPM0'],4)
xlabel('Chi-squared','FontSize',12)
ylabel('Freuquency','FontSize',12)
title('Distributions','FontSize',16)

ne   = sum(SPM > max(SPM0));
pe   = 1 - spm_Icdf(ne,length(e),1/length(e));
fprintf('\nOminbus p-value: p < %-0.5f\n',pe);

return
 

% NOTES: Phase portrait and entropies (temperature)
%==========================================================================
spm_figure('GetWin','Phases and entropies');
T     = 1024;
t     = (T - 512):T;
n     = 16;
Hp    = zeros(1,n);
Hq    = zeros(1,n);
Pe    = linspace(-4,4,n);
Pp    = P;
for i = 1:n
    
    
    % parameters (energy and state-dependent forces)
    %------------------------------------------------------------------
    Pp.e    = exp(Pe(i));
    
    % entropy of dynamics
    %------------------------------------------------------------------
    [Q,X]   = spm_Manifold_solve(x,u,Pp,T,1/32,0);
    
    
    % entropy of dynamics
    %------------------------------------------------------------------
    q       = squeeze(Q(1,:,t));
    Hq(i) = spm_logdet(cov(q'));
    
    % configuration entropy
    %------------------------------------------------------------------
    p       = squeeze(X(1,:,t));
    Hp(i) = spm_logdet(cov(p'));
    
    % update with graphics
    %------------------------------------------------------------------
    subplot(2,2,1)
    plot(Pe,Hp,Pe,Hq)
    title('Phases and entropy','FontSize',16)
    xlabel('Entropy','FontSize',12)
    ylabel('Potential energy','FontSize',12)
    legend({'Structural','Functional'})
    axis square
    
    subplot(2,2,2)
    plot(Pe,Hp - Hq)
    title('Difference','FontSize',16)
    xlabel('Entropy','FontSize',12)
    ylabel('Potential energy','FontSize',12)
    axis square
    
    % plot dynamics
    %--------------------------------------------------------------------------
    subplot(2,2,3)
    plot(1:T,squeeze(X(1,:,:)),':',t,squeeze(X(1,:,t)),'-')
    axis([1 T -64 64])
    xlabel('time','FontSize',12)
    title('Dynamics','FontSize',16)
    axis square
    
    % plot entropies
    %--------------------------------------------------------------------------
    subplot(2,2,4)
    plot(Hq(:),Hp(:),'.',Hq(:),Hp(:),'o')
    xlabel('Functional entropy','FontSize',12)
    ylabel('Structural entropy','FontSize',12)
    title('Entropies','FontSize',16)
    axis square
    
end

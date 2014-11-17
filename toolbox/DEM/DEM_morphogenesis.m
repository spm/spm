function DEM = DEM_morphogenesis
% This demonstration uses active inference (as implemented in spm_ALAP) to
% illustrate birdsong and communication using predictive coding. In this
% example, priors on high-level sensorimotor constructs (e.g., in the
% avian higher vocal centre) are used to generate proprioceptive
% predictions (i.e., motor commands) so that the bird can sing. However, in
% the absence of sensory attenuation, the slight differences between
% descending predictions and the sensory consequences of self-made songs
% confound  singing. This means that sensory attenuation is required so
% that the bird can either sing or listen.  By introducing a second bird
% and alternating between singing and listening respectively, one can
% simulate communication through birdsong. This is illustrated with one
% cycle of singing and listening, where the high level expectations about
% hidden states become synchronised; in effect, the two birds are singing
% from the same 'hymn sheet' or narrative and can be regarded as
% communicating in the sense of pragmatics. The first bird's expectations
% are shown in red, while the second bird's are shown in cyan.
%
% To simulate learning of each other's (high-level) attractor, set LEARN to
% one in the  main script.. To separate the birds – and preclude
% communication (or synchronisation chaos) set NULL to 1.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_morphogenesis.m 6263 2014-11-17 13:48:36Z karl $
 

% preliminaries
%--------------------------------------------------------------------------
rng('default')

N        = 64;                             % length of process (bins)

% generative process and model
%==========================================================================
M(1).E.d = 1;                              % approximation order
M(1).E.n = 2;                              % embedding order
M(1).E.s = 1;                              % smoothness

% priors (prototype)
%--------------------------------------------------------------------------
p(:,:,1) = [1 1 1 1];
p(:,:,2) = [1 1 0 0];
p(:,:,3) = [0 1 0 0];

[y x] = find(p(:,:,1));
P.x   = spm_detrend([x' y'])'; 

% signalling of each cell type
%--------------------------------------------------------------------------
for i = 1:size(p,3);
    P.s(i,:) = spm_vec(p(:,:,i))';
end
P.s   = double(P.s);
P.c   = morphogenesis(P.x,P.s);           % signal sensed at each position

% initialise action and expectations
%--------------------------------------------------------------------------
n     = size(P.x,2);
v.i   = randn(n,n)/8;                     % states (identity)
v.q   = zeros(3,n);                       % states (expression)

g     = Mg([],v,P);
a.x   = randn(size(P.x))/8;               % action (chemotaxis)
a.s   = g.s;                              % action (signal release)



% generative process 
%==========================================================================

% level 1 of generative process
%--------------------------------------------------------------------------
G(1).g  = @(x,v,a,P) Gg(x,v,a,P);
G(1).V  = exp(16);                         % precision (noise)
G(1).U  = exp(2);                          % precision (action)
G(1).pE = a;                               % form (action)

 
% level 2; causes (action)
%--------------------------------------------------------------------------
G(2).a  = spm_vec(a);                      % endogenous cause (action)
G(2).v  = 0;                               % exogenous  cause
G(2).V  = exp(16);


% generative model
%==========================================================================

% level 1 of the generative model: 
%--------------------------------------------------------------------------
M(1).g  = @(x,v,P) Mg([],v,P);
M(1).W  = exp(14);
M(1).v  = g;
M(1).V  = exp(4);
M(1).pE = P;

% level 2: 
%--------------------------------------------------------------------------
M(2).v  = v;
M(2).V  = diag([(zeros(n*n,1) + exp(-4)); (zeros(3*n,1) + exp(4))]);


% hidden cause and prior expectations
%--------------------------------------------------------------------------
for t = 1:N
    for i = 1:3
        q(i,:) = zeros(n,1) + spm_phi(t - (i - 1)*N/4);
    end
    u.i    = zeros(n,n);                   % states (identity)
    u.q    = q;                            % states (expression)
    U(:,t) = spm_vec(u);
end
C     = zeros(1,N);

% assemble model structure
%--------------------------------------------------------------------------
DEM.M = M;
DEM.G = G;
DEM.C = C;
DEM.U = U;

% reset initial hidden states and invert
%==========================================================================
DEM   = spm_ADEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU)


% graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
subplot(2,2,3); cla;
A     = max(abs(P.x(:)))*3/2;

for t = 1:N
    v = spm_unvec(DEM.qU.a{2}(:,t),a);
    for i = 1:n
        x = v.x(1,i);
        c = v.s(end - 2:end,i);
        c = max(min(c,1),0);
        plot(t,x,'.','markersize',t + 4,'color',c); hold on
    end
end

% target morphology
%--------------------------------------------------------------------------
for i = 1:n
    x = P.x(1,i);
    c = P.s(end - 2:end,i);
    c = max(min(c,1),0);
    plot(N + 4,x,'.','markersize',t + 4,'color',c); hold on
end

title('morphogenesis','Fontsize',16)
xlabel('time')
ylabel('llocation')
set(gca,'Color','k');
axis square, box off

% expected signal concentrations
%--------------------------------------------------------------------------
y     = linspace(min(P.x(1,:)),max(P.x(1,:)),32);
for i = 1:size(y,2)
    c = morphogenesis(P.x,P.s,y(:,i));
    c = c(end - 2:end);
    c = max(min(c,1),0);
    plot(N + 8,y(:,i),'.','markersize',t + 8,'color',c); hold on
end

title('morphogenesis','Fontsize',16)
xlabel('time')
ylabel('llocation')
set(gca,'Color','k','YLim',[-1 1]*A);
axis square, box off


% movies
%--------------------------------------------------------------------------
subplot(2,2,1); cla;
for t = 1:N
    v = spm_unvec(DEM.qU.a{2}(:,t),a);
    hold off
    for i = 1:n
        x = v.x(:,i);
        c = v.s(end - 2:end,i);
        c = max(min(c,1),0);
        plot(x(1),x(2),'.','markersize',t + 4,'color',c); hold on
    end
    set(gca,'Color','k');
    axis square, box off
    axis([-1 1 -1 1]*A)
    drawnow
    
    % save
    %------------------------------------------------------------------
    Mov(t) = getframe(gca);
    
end

set(gca,'Userdata',{Mov,8})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Extrinsic (left click for movie)','FontSize',16)
xlabel('location')

% target morphology
%--------------------------------------------------------------------------
subplot(2,2,2); cla
for i = 1:n
    x = P.x(:,i);
    c = P.s(end - 2:end,i);
    c = max(min(c,1),0);
    plot(x(1),x(2),'.','markersize',t + 4,'color',c); hold on
end

title('morphogenesis','Fontsize',16)
xlabel('time')
ylabel('llocation')
set(gca,'Color','k');
axis([-1 1 -1 1]*A)
axis square, box off
hold off
    

return


% Equations of motion and observer functions
%==========================================================================

% sensed signal
%--------------------------------------------------------------------------
function c = morphogenesis(x,s,y)
% x - location of cells
% s - signals released
% y - location of sampling [default: x]
%__________________________________________________________________________

% preliminaries
%--------------------------------------------------------------------------
if nargin < 3; y = x; end                  % sample locations
n     = size(y,2);                         % number of locations
m     = size(s,1);                         % number of signals
k     = [1 1 1 1]';                        % signal decay over space 
c     = zeros(m,n);                        % signal sensed at each location
k     = k(1:m);
K     = [1 0;0 0];
for i = 1:n
    for j = 1:size(x,2)
        
        % distance
        %------------------------------------------------------------------
        d      = y(:,i) - x(:,j);
        d      = sqrt(d'*K*d);
        
        % signal concentration
        %------------------------------------------------------------------
        c(:,i) = c(:,i) + exp(-k*d).*s(:,j);

    end
end


% first level process: generating input
%--------------------------------------------------------------------------
function g = Gg(x,v,a,P)

a     = spm_unvec(a,P);
g.x   = a.x;                             % position
g.s   = a.s;                             % intrinsic signal
g.c   = morphogenesis(a.x,a.s);          % extrinsic signal

% first level model: mapping hidden causes to sensations
%--------------------------------------------------------------------------
function g = Mg(x,v,P)

p    = spm_softmax(v.i);                 % expected indentity
g.x  =  P.x*p;                           % position
g.s  = (P.s*p).*v.q;                     % intrinsic signal
g.c  = (P.c*p).*v.q;                     % extrinsic signal




 
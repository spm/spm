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
% $Id: DEM_morphogenesis.m 6254 2014-11-04 18:24:21Z karl $
 

% preliminaries
%--------------------------------------------------------------------------
rng('default')

n   = 5;                                   % number of cells
m   = 3;                                   % number of signals
N   = 64;                                  % length of stimulus (bins)

% generative process and model
%==========================================================================
M(1).E.d = 1;                              % approximation order
M(1).E.n = 2;                              % embedding order
M(1).E.s = 1/2;                            % smoothness

% priors (prototype)
%--------------------------------------------------------------------------
P.x     = (1:n)/2;                         % position of each cell
P.s     = [(P.x > 0);                      % signalling of each cell
           (P.x > mean(P.x));
           (P.x < mean(P.x));
           (P.x == mean(P.x))];
P.x     = spm_detrend(P.x')';
P.s     = double(P.s);
P.c     = morphogenesis(P.x,P.s);          % signal sensed at each position

% initialise action and expectations
%--------------------------------------------------------------------------
for i = 1:n
    x(i).i = randn(n,1)/2;                 % states (identity)
end
g       = Mg(x,[],P);
a.x     = zeros(1,n);                      % states (position)
a.s     = g.s;                             % states (chemotaxis)



% generative process 
%==========================================================================

% level 1 of generative process
%--------------------------------------------------------------------------
G(1).g  = @(x,v,a,P) Gg(x,v,a,P);
G(1).V  = exp(16);                         % precision (noise)
G(1).U  = exp(4);                          % precision (action)
G(1).pE = a;                               % precision (action)

 
% level 2; causes (action)
%--------------------------------------------------------------------------
G(2).a  = spm_vec(a);                      % endogenous cause (action)
G(2).v  = 0;                               % exogenous  cause
G(2).V  = exp(16);


% generative model
%==========================================================================

% level 1 of the generative model: 
%--------------------------------------------------------------------------
M(1).f  = @(x,v,P) Mf(x,v,P);
M(1).g  = @(x,v,P) Mg(x,v,P);
M(1).x  = x;
M(1).W  = exp(8);
M(1).V  = exp(4);
M(1).pE = P;

% level 2: 
%--------------------------------------------------------------------------
M(2).v  = 0;
M(2).V  = exp(16);


% hidden cause and prior expectations
%--------------------------------------------------------------------------
C     = zeros(1,N);

% assemble model structure
%--------------------------------------------------------------------------
DEM.M = M;
DEM.G = G;
DEM.C = C;
DEM.U = C;

% reset initial hidden states and invert
%==========================================================================
DEM   = spm_ADEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU)


% graphics
%--------------------------------------------------------------------------
subplot(2,2,3); cla;

for t = 1:N
    v = spm_unvec(DEM.qU.a{2}(:,t),a);
    for i = 1:n
        x = v.x(:,i);
        c = v.s(end - 2:end,i);
        c = max(min(c,1),0);
        plot(t,x,'.','markersize',t + 4,'color',c); hold on
    end
end

% target morphology
%--------------------------------------------------------------------------
for i = 1:n
    x = P.x(:,i);
    c = P.s(end - 2:end,i);
    c = max(min(c,1),0);
    plot(N + 4,x,'.','markersize',t + 4,'color',c); hold on
end

title('morphogenesis','Fontsize',16)
xlabel('time')
ylabel('llocation')
set(gca,'Color','k');
axis square, box off

return

% illustrate responses with sonogram (and sound file)
%==========================================================================
% spm_figure('GetWin','Figure 1'); clf




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
k     = [4 2 1 1]';                         % signal decay over space 
c     = zeros(size(s,1),size(y,2));        % signal sensed at each location
for i = 1:size(y,2)
    for j = 1:size(x,2)
        
        % distance
        %------------------------------------------------------------------
        d      = y(:,i) - x(:,j);
        d      = sqrt(d'*d);
        
        % signal concentration
        %------------------------------------------------------------------
        c(:,i) = c(:,i) + exp(-k*abs(d)).*s(:,j);

    end
end


% first level process: generating input
%--------------------------------------------------------------------------
function g = Gg(x,v,a,P)
a     = spm_unvec(a,P);
g.s   = a.s;
g.c   = morphogenesis(a.x,a.s);

% first level model: mapping hidden causes to sensations
%--------------------------------------------------------------------------
function g = Mg(x,v,P)
n     = length(x);
for i = 1:n 
    p        = spm_softmax(x(i).i);
    g.s(:,i) = P.s*p;
    g.c(:,i) = P.c*p;
end

% first level model: flow of hidden states
%--------------------------------------------------------------------------
function f = Mf(x,v,P)
n   = length(x);
f   = x;
for i = 1:n
    f(i).i = ((speye(n,n) - 1)*( 1./(1 + exp(-x(i).i))) - x(i).i/8  + 1)/4;
end


 
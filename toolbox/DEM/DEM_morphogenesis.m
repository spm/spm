function DEM_morphogenesis
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
% $Id: DEM_morphogenesis.m 6235 2014-10-12 10:03:05Z karl $
 

% preliminaries
%--------------------------------------------------------------------------
rng('default')

A   = 4;                                   % number of cells
N   = 64;                                 % length of stimulus (bins)

% generative process and model
%==========================================================================
M(1).E.d        = 1;                       % approximation order
M(1).E.n        = 2;                       % embedding order
M(1).E.s        = 1/2;                     % smoothness
 

% initialise parameters, states and precisions
%--------------------------------------------------------------------------
x     = cell(A,1);
a     = cell(A,1);

% level 1 of generative process
%--------------------------------------------------------------------------
for i = 1:A
    x{i}.x = randn(1,1)/4;                 % states (position)
    x{i}.v = zeros(1,1)/4;                 % states (velocity)
    x{i}.c = spm_softmax(randn(2,1));      % states (chemotaxis)
end
G(1).f  = @(x,v,a,P) Gf(x,v,a,P);
G(1).g  = @(x,v,a,P) Gg(x,v,a,P);          % 
G(1).x  = x;                               % hidden state
G(1).V  = exp(16);                         % precision (noise)
G(1).W  = exp(16);                         % precision (states)
G(1).U  = exp(4);                          % precision (action)

 
% level 2; causes
%--------------------------------------------------------------------------
for i = 1:A
    a{i}.x = zeros(1,1);                   % action (migration)
    a{i}.c = zeros(2,1);                   % action (release)
end
G(2).v  = 0;                               % exogenous  cause
G(2).a  = spm_vec(a);                      % endogenous cause (action)
G(2).V  = exp(16);



% level 1 of the generative model: 
%--------------------------------------------------------------------------
clear x
for i = 1:A
    x{i}.i = zeros(3,1);                   % states (identity)
end
M(1).f  = @(x,v,P) Mf(x,v,P);
M(1).g  = @(x,v,P) Mg(x,v,P);
M(1).x  = x;
M(1).W  = exp(4);
M(1).V  = exp(4);
M(1).pE = x; 

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
    
return

% illustrate responses with sonogram (and sound file)
%==========================================================================
% spm_figure('GetWin','Figure 1'); clf
% 
% 
% for t = 1:T
%     x = DEM.qU.x{2}([1 2 3],:);
%     y = DEM{t}.qU.x{2}([1 2 3] + 3,:);
% end
% 
% plot(x',y')
% title('Synchronization manifold','Fontsize',16)
% xlabel('second level expectations (first bird)')
% ylabel('second level expectations (second bird)')
% axis square



% Equations of motion and observer functions
%==========================================================================


% first level process: Genersting input
%--------------------------------------------------------------------------
function g = Gg(x,v,a,P)
n     = length(x);
m     = length(x{1}.c);
for i = 1:n
    g(i).c = 0;
    g(i).s = 0;
    for j = 1:n
        
        % distance
        %------------------------------------------------------------------
        d      = x{i}.x - x{j}.x;
        
        % concentration (local)
        %------------------------------------------------------------------
        q      = exp(-8*abs(d));
        g(i).c = g(i).c + q;
        
        % concentration (global)
        %------------------------------------------------------------------
        q      = exp(-1*abs(d));
        for k = 1:m
           c(k,1) = x{j}.c(k)*q;
        end
        g(i).s = g(i).s + c;

    end
end

% first level model: mapping hidden causes to sensations
%--------------------------------------------------------------------------
function g = Mg(x,v,P)
n     = length(x);
m     = 2;
for i = 1:n 
    g(i).c = 1 + 1/2;
    g(i).s = zeros(2,1);
end

% first level process: Hidden states
%--------------------------------------------------------------------------
function f = Gf(x,v,a,P)
n     = length(x);
m     = length(x{1}.c);
for i = 1:n
    q{i}.x = x{i}.x;
    q{i}.c = x{i}.c;
end
a   = spm_unvec(a,q);
f   = x;
for i = 1:n
    f{i}.x = x{i}.v;
    f{i}.v = a{i}.x/16 - 4*x{i}.v;
    f{i}.c = a{i}.c/16 - x{i}.c/256;
end


% first level model: flow of hidden states
%--------------------------------------------------------------------------
function f = Mf(x,v,P)
n   = length(x);
f   = x;
for i = 1:n
    f{i}.i = 0*x{i}.i;
end

 
function ADEM_learning
% Value learning demo using the mountain car problem. This demo questions
% the need for reinforcement learning and related paradigms from
% machine-learning, when trying to optimise the behaviour of an agent.  We
% show that it is fairly simple to teach an agent complicated and adaptive
% behaviours under the free-energy principle.  This principle suggests that
% agents adjust their internal states and sampling of the environment to
% minimize their free-energy.  In this context, free-energy represents a
% bound on the probability of being in a particular state, given the nature
% of the agent, or more specifically the model of the environment an agent
% entails.  We show that such agents learn causal structure in the
% environment and sample it in an adaptive and self-supervised fashion.
% The result is a behavioural policy that reproduces exactly the policies
% that are optimised by reinforcement learning and dynamic programming.
% Critically, at no point do we need to invoke the notion of reward, value
% or utility.

%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_learning.m 6801 2016-05-29 19:18:06Z karl $
 
 
% generative model
%==========================================================================
rng('default')
DEMO     = 0;                           % switch for demo

G(1).E.s = 1/2;                         % smoothness
G(1).E.n = 4;                           % smoothness
G(1).E.d = 2;                           % smoothness
 
% parameters
%--------------------------------------------------------------------------
P.a     = 0;
P.b     = [0 0];
P.c     = [0 0 0 0];
P.d     = 0;
P0      = P;
pC      = speye(length(spm_vec(P)));
pC(end) = 0;
 
% level 1
%--------------------------------------------------------------------------
G(1).x  = [0; 0];
G(1).f  = 'spm_fx_mountaincar';
G(1).g  = @(x,v,a,P) x;
G(1).pE = P;
G(1).pC = pC;
G(1).V  = exp(8);                       % error precision
G(1).W  = exp(8);                       % error precision
 
% level 2
%--------------------------------------------------------------------------
G(2).a  = 0;                            % action
G(2).v  = 0;                            % inputs
G(2).V  = exp(16);
G       = spm_ADEM_M_set(G);
 
 
% desired equilibrium density and state space
%==========================================================================
G(1).fq = 'spm_mountaincar_Q';
 
% create X - coordinates of evaluation grid
%--------------------------------------------------------------------------
nx      = 32;
x{1}    = linspace(-1,1,nx);
x{2}    = linspace(-1,1,nx);
[X,x]   = spm_ndgrid(x);
G(1).X  = X;
 
 
% optimise parameters so that p(y|G) maximises the cost function
%==========================================================================
 
% optimise parameters: (NB an alternative is P = spm_fp_fmin(G));
%--------------------------------------------------------------------------
if DEMO
    
    % Optimise parameters of fictive forces using KL control
    %----------------------------------------------------------------------
    P.a = 2.8985;
    P.b = [-1.7468 3.3942];
    P.c = [-1.0146 -0.0021 -4.7885 -0.6077];
    P.d = 0;
    
    if DEMO > 1
        P   = spm_fmin('spm_mountaincar_fun',P,pC,G);
        P.d = 0;
    end
    
    G(1).pE = P;
    disp(P)
    save mountaincar_model G
end
 
 
% or load previously optimised environment
%--------------------------------------------------------------------------
load mountaincar_model
P     = G(1).pE;
 
% plot flow fields and nullclines
%==========================================================================
spm_figure('GetWin','Figure 1');
 
nx    = 64;
x{1}  = linspace(-2,2,nx);
x{2}  = linspace(-2,2,nx);
M     = G;
 
% uncontrolled flow (P0)
%--------------------------------------------------------------------------
M(1).pE = P0;
subplot(3,2,1)
spm_fp_display_density(M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('flow and equilibrium density','Fontsize',16)
 
subplot(3,2,2)
spm_fp_display_nullclines(M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('nullclines','Fontsize',16)
 
% controlled flow (P0)
%--------------------------------------------------------------------------
M(1).pE = P;
subplot(3,2,3)
spm_fp_display_density(M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('controlled','Fontsize',16)
 
subplot(3,2,4)
spm_fp_display_nullclines(M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('controlled','Fontsize',16)
drawnow
 

% recognition model: learn the controlled environmental dynamics
%==========================================================================
M       = G;
M(1).g  = @(x,v,P)x;
 
% make a niave model (M)
%--------------------------------------------------------------------------
M(1).pE = P0;
M(1).pC = exp(8);
M(1).V  = exp(8);
M(1).W  = exp(8);

% teach naive model by exposing it to a controlled environment (G)
%--------------------------------------------------------------------------
clear DEM
 
% perturbations
%--------------------------------------------------------------------------
n     = 16;
i     = (1:n)*32;
C     = sparse(1,i,randn(1,n));
C     = spm_conv(C,4);
 
DEM.M = M;
DEM.G = G;
DEM.C = C;
DEM.U = C;
 
% optimise recognition model
%--------------------------------------------------------------------------
if DEMO
    DEM.M(1).E.nE = 16;
    DEM           = spm_ADEM(DEM);
    save mountaincar_model G DEM
end

load mountaincar_model

 
% replace priors with learned conditional expectation
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
 
M(1).pE = DEM.qP.P{1};
M(1).pC = [];

subplot(3,2,5)
spm_fp_display_density(M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('learnt','Fontsize',16)
 
subplot(3,2,6)
spm_fp_display_nullclines(M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('learnt','Fontsize',16)

 
% evaluate performance under active inference
%==========================================================================
 
% create uncontrolled environment (with action)
%--------------------------------------------------------------------------
G(1).pE   = P0;
G(1).pE.d = 1;
G(1).U    = exp(8);
G(1).V    = exp(16);
G(1).W    = exp(16);

% create DEM structure (and perturb the real car with fluctuations)
%--------------------------------------------------------------------------
N       = 128;
U       = sparse(1,N);
C       = spm_conv(randn(1,N),8)/4;      % pertubations
DEM.G   = G;
DEM.M   = M;
DEM.C   = U;
DEM.U   = U;
DEM     = spm_ADEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM');
spm_DEM_qU(DEM.qU,DEM.pU)
 
subplot(2,2,3)
spm_fp_display_nullclines(M,x);hold on
plot(DEM.pU.v{1}(1,:),DEM.pU.v{1}(2,:),'b'), hold off
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('learnt','Fontsize',16)

% movie
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
clf, subplot(3,1,2)
drawnow
spm_mountaincar_movie(DEM)
title('click car for movie','FontSize',16)




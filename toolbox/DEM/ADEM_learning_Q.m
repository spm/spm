% Value learning using the mountain car problem: This demo explores the
% use of the Helmholtz (Ao) decomposition directly as a generative model
% for flow. The advantage of this is that the value function (log ergodic
% density) can be specified directly in terms of the desired density. This
% means the divergence-free and diffusion parameters can then be learned to
% accommodate environmental constraints on flow. The constraints on value
% are that is has to have maximum at and only at) the desired location.
%
% This scheme was not pursued to closure
 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_learning_Q.m 4516 2011-10-07 19:18:32Z karl $
 
 
% generative model
%==========================================================================
DEMO     = 1;                          % switch for demo
G(1).E.s = 1/2;                        % smoothness
G(1).E.n = 6;                          % embedding dimension
G(1).E.d = 2;                          % embedding dimension
 
% parameters
%--------------------------------------------------------------------------
P.V     = zeros(6,1);
P.Q     = zeros(6,1);
P.J     = [-8 -2];
pC      = speye(length(spm_vec(P)));
 
% level 1
%--------------------------------------------------------------------------
G(1).x  = [0; 0];
G(1).f  = 'spm_fx_mountaincar';
G(1).g  = inline('x','x','v','a','P');
G(1).V  = exp(16);                           % error precision
G(1).W  = diag([exp(16) exp(8)]);            % error precision
 
% level 2
%--------------------------------------------------------------------------
G(2).a  = 0;                                % action
G(2).v  = 0;                                % inputs
G(2).V  = exp(16);
G       = spm_ADEM_M_set(G);
 
 
% recognition model
%==========================================================================
M       = G;
M(1).f  = 'spm_fx_mountaincar_Q';
M(1).g  = inline('x','x','v','P');
M(1).pE = P;
M(1).pC = exp(0);
 
% examine the true and prior ergodic densities under these flows
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
nx      = 32;
x{1}    = linspace(-2,2,nx);
x{2}    = linspace(-2,2,nx);
[X x]   = spm_ndgrid(x);
 
subplot(2,2,1)
FG   = spm_fp_display_density(G,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title( 'true',  'Fontsize',16)
 
subplot(2,2,2)
FM   = spm_fp_display_density(M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title( 'prior policy',  'Fontsize',16)
 
subplot(2,2,3)
plot(FG)
xlabel('position','Fontsize',12)
ylabel('flow','Fontsize',12)
title( 'true',  'Fontsize',16)
axis square
 
subplot(2,2,4)
plot(FM)
xlabel('position','Fontsize',12)
ylabel('flow','Fontsize',12)
title( 'prior policy',  'Fontsize',16)
axis square
drawnow
 
% known (deterministic) exploratory perturbations
%--------------------------------------------------------------------------
clear DEM
N     = 512;
DEM.M = M;
DEM.G = G;
DEM.C = spm_conv(randn(1,N),16);
DEM.U = DEM.C;
 
% optimise recognition model
%--------------------------------------------------------------------------
DEM.M(1).E.nE = 4;
if DEMO
    
    for i = 1:4
        
        % initialise position
        %------------------------------------------------------------------
        x0         = [2*rand - 1; 0];
        DEM.M(1).x = x0;
        DEM.G(1).x = x0;
        DEM.G(2).a = 0;
        DEM        = spm_ADEM(DEM);
        DEM        = spm_ADEM_update(DEM);
        
        spm_figure('GetWin','Figure 2');
        spm_DEM_qU(DEM.qU,DEM.pU)
        
    end
    
    % save
    %----------------------------------------------------------------------
    save DEM_Q DEM
    
else
    
    % load
    %----------------------------------------------------------------------
    load DEM_Q
end
 
 
% empirical priors (policy)
%==========================================================================
 
% create X - coordinates of evaluation grid
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');clf
nx      = 64;
x{1}    = linspace(-2,2,nx);
x{2}    = linspace(-2,2,nx);
[X x]   = spm_ndgrid(x);
 
DEM.M(1).pE = DEM.qP.P{1};
DEM.M(1).pC = [];
 
subplot(2,1,1)
spm_fp_display_density(DEM.M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title( 'policy',  'Fontsize',16)
 
 
% evaluate performance under active inference
%==========================================================================
 
% create DEM structure
%--------------------------------------------------------------------------
N       = 128;
U       = sparse(1,N);
DEM.C   = U;
DEM.U   = U;
 
x0         = [0; -1/2];
DEM.M(1).x = x0;
DEM.G(1).x = x0;
DEM.G(2).a = 0;
DEM        = spm_ADEM(DEM);
 
% plot solutions
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_DEM_qU(DEM.qU,DEM.pU)
 
subplot(2,2,3)
plot(DEM.pU(1).x{1}(1,:),DEM.pU(1).x{1}(2,:))
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('solutions','Fontsize',16)
axis square
drawnow
 
for i = 1:16
    
    % initialise position
    %----------------------------------------------------------------------
    x0          = [rand - 1/2; rand - 1/2];
    DEM.M(1).x = x0;
    DEM.G(1).x = x0;
    DEM.G(2).a = 0;
    DEM        = spm_ADEM(DEM);
    
    % overlay true values
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 3');
    subplot(2,1,2)
    hold on
    plot(DEM.pU(1).x{1}(1,:),DEM.pU(1).x{1}(2,:))
    xlabel('position','Fontsize',12)
    ylabel('velocity','Fontsize',12)
    title('solutions','Fontsize',16)
    axis square
    drawnow
    
end
 
 
 



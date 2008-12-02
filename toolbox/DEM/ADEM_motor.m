% This demo illustrates how action can fulfil prior expectations by
% explaining away sensory prediction errors prescribed by desired movement
% trajectory. It is based on the same linear convolution model of the
% motor plant considered in the visual tracking example. Here, we induce
% prediction errors; not through exogenous perturbation to sensory input
% but through tight priors encoding a desired or expected trajectory. We 
% then show how the movement is robust to changes in the true motor
% dynamics and other exogenous perturbations, late in movement execution
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_motor.m 2521 2008-12-02 19:49:39Z karl $
 
% Recognition model (linear for expediency)
%==========================================================================
clear
M(1).E.linear = 1;                          % linear model
M(1).E.s      = 1/2;                        % smoothness
M(1).E.n      = 6;                          % smoothness
M(1).E.d      = 2;                          % smoothness
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
pE.f    = [-1  4   ;                        % the Jacobian for the
           -2 -1]/8;                        % hidden sates
pE.g    = [spm_dctmtx(4,2)]/8;              % the mixing parameters
pE.h    = [1; 0];                           % input parameter
M(1).x  = [0; 0];
M(1).f  = inline('P.f*x + P.h*v','x','v','P');
M(1).g  = inline('P.g*x','x','v','P');
M(1).pE = pE;                               % prior expectation
M(1).V  = exp(4);                           % error precision
M(1).W  = exp(8);                           % error precision
 
% level 2: with informative priors on movement
%--------------------------------------------------------------------------
M(2).v  = 0;                                % inputs
M(2).V  = exp(16);
 
% generative model
%==========================================================================
G       = M;
pE.a    = [1; 0];                           % action parameter
 
% first level
%--------------------------------------------------------------------------
G(1).f  = inline('P.f*x + P.h*v + P.a*a','x','v','a','P');
G(1).g  = inline('P.g*x','x','v','a','P');
G(1).pE = pE;                               % prior expectation
G(1).V  = exp(16);                           % error precision
G(1).W  = exp(16);                           % error precision
 
% second level
%--------------------------------------------------------------------------
G(2).a  = 0;                                % action
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N       = 64;                                 % length of data sequence
C       = exp(-([1:N] - 12).^2/(4.^2));       % this is the prior cause;
DEM.G   = G;
DEM.M   = M;
DEM.C   = sparse(1,N);
DEM.U   = C;
DEM0    = spm_ADEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM0.qU,DEM0.pU)
 
 
% repeat with a late perturbation
%--------------------------------------------------------------------------
DEM1    = DEM;
DEM1.C  = -exp(-([1:N] - 18).^2/(2.^2))/2;       % this is the prior cause;
DEM1    = spm_ADEM(DEM1);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM1.qU,DEM1.pU)
 
 
 
% repeat with twice the motor gain (P.a)
%--------------------------------------------------------------------------
DEM2            = DEM;
DEM2.G(1).pE.a  = pE.a*2;
DEM2            = spm_ADEM(DEM2);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM2.qU,DEM2.pU)
 
 
% plot results
%==========================================================================
spm_DEM_qU(DEM0.qU,DEM0.pU)
subplot(2,2,2)
title('desired occulomotor state','FontSize',16)
 
 
spm_figure('GetWin','Graphics');
 
% canonical
%--------------------------------------------------------------------------
subplot(3,2,1)
plot(DEM0.pU.x{1}(1,:),DEM0.pU.x{1}(2,:),DEM0.qU.x{1}(1,:),DEM0.qU.x{1}(2,:),':')
axis square
title('desired trajectory','FontSize',16)
xlabel('displacement','FontSize',14)
legend({'real', 'perceived'})
axis([-2 2 -2 2])
 
subplot(3,2,2)
plot([1:N],DEM0.qU.a{2},[1:N],DEM0.qU.v{2},':',[1:N],DEM0.pU.v{2},'-.')
axis square
title('action and causes','FontSize',16)
xlabel('displacement','FontSize',14)
legend({'action', 'perceived cause', 'exogenous cause'})
axis([1 N -1/2 1])
 
% late perturbation
%--------------------------------------------------------------------------
subplot(3,2,3)
plot(DEM1.pU.x{1}(1,:),DEM1.pU.x{1}(2,:),DEM1.qU.x{1}(1,:),DEM1.qU.x{1}(2,:),':')
axis square
title('with perturbation','FontSize',16)
xlabel('displacement','FontSize',14)
axis([-2 2 -2 2])
 
subplot(3,2,4)
plot([1:N],DEM1.qU.a{2},[1:N],DEM1.qU.v{2},':',[1:N],DEM1.pU.v{2},'-.')
axis square
title('action and causes','FontSize',16)
xlabel('time','FontSize',14)
axis([1 N -1/2 1])
 
% change in motor gain
%--------------------------------------------------------------------------
subplot(3,2,5)
plot(DEM2.pU.x{1}(1,:),DEM2.pU.x{1}(2,:),DEM2.qU.x{1}(1,:),DEM2.qU.x{1}(2,:),':')
axis square
title('change in motor gain','FontSize',16)
xlabel('displacement','FontSize',14)
axis([-2 2 -2 2])
 
subplot(3,2,6)
plot([1:N],DEM2.qU.a{2},[1:N],DEM2.qU.v{2},':',[1:N],DEM2.pU.v{2},'-.')
axis square
title('action and causes','FontSize',16)
xlabel('time','FontSize',14)
axis([1 N -1/2 1])

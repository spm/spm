% DEM demo for active inference (i.e. action-perception optimisation of free
% energy).  This simulation calls on spm_ADEM to simulate visual sampling of
% the world and demonstrate retinal stabilisation or visual tracking. We
% simulate a simple 2-D plaid stimulus and subject it to an exogenous
% perturbations. By employing tight and broad priors on the location of the
% stimulus, we can show that action does and does not explain away the visual
% consequences of the perturbation (i.e., the movement is seen or not).  This
% illustrates how one can reframe stabilisation or tracking in terms of
% sampling sensory input to ensure conditional expectations are met; and
% how these expectations are shaped by prior expectations.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_visual.m 3695 2010-01-22 14:18:14Z karl $
 
 
% recognition model (M)
%==========================================================================
clear
M(1).E.s      = 1/2;                        % smoothness
M(1).E.n      = 6;                          % smoothness
M(1).E.d      = 2;                          % smoothness
 
% level 1:
% the hidden states correspond to horizontal and vertical displacement of a
% Gaussian modulated plaid
%--------------------------------------------------------------------------
pE.f    = [-1  4   ;                        % the Jacobian for the
           -2 -1]/8;                        % hidden states (x and y)
pE.h    = [1; 0];                           % input parameters
M(1).x  = [0; 0];
M(1).f  = inline('P.f*x + P.h*v','x','v','P');
M(1).g  = inline('ADEM_plaid(x)','x','v','P');
M(1).pE = pE;                               % prior expectation
M(1).V  = exp(4);                           % error precision
M(1).W  = exp(8);                           % error precision
 
% level 2:
% a perturbation to x(1)
%--------------------------------------------------------------------------
M(2).v  = 0;                                 % inputs
M(2).V  = exp(-16);                          % flat priors on movement
 
% Generative model (G)
%==========================================================================
G       = M;
pE.a    = [1; 0];                           % action parameter
 
% first level
%--------------------------------------------------------------------------
G(1).f  = inline('P.f*x + P.h*v + P.a*a','x','v','a','P');
G(1).g  = inline('ADEM_plaid(x)','x','v','a','P');
G(1).pE = pE;                               % prior expectation
G(1).V  = exp(16);                          % error precision
G(1).W  = exp(8);                           % error precision
 
% second level
%--------------------------------------------------------------------------
G(2).a  = 0;                                % action
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N       = 64;                               % length of data sequence
C       = exp(-([1:N] - 12).^2/(4.^2));     % this is the perturbation
DEM.G   = G;
DEM.M   = M;
DEM.C   = C;
DEM     = spm_ADEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
 
 
% repeat with informative priors
%--------------------------------------------------------------------------
DEM.M(1).V  = exp(4);
DEM.M(1).W  = exp(8);
DEM.M(2).V  = exp(16);
ADEM        = spm_ADEM(DEM);
spm_DEM_qU(ADEM.qU,ADEM.pU)
 
% plot results
%==========================================================================
spm_figure('GetWin','Graphics');
x   = [-4 4];
 
subplot(3,1,1)
imagesc(x,x,ADEM_plaid([0; 0],128))
axis square
title('stimulus','FontSize',16)
x   = linspace(-3.5,3.5,6);
hold on
plot(kron(x.^0,x),kron(x,x.^0),'.w','MarkerSize',16), hold off
 
subplot(3,2,3)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:),DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),':')
axis square
title('under flat priors','FontSize',16)
xlabel('displacement','FontSize',14)
legend({'real', 'perceived'})
axis([-1 1 -1 1]*2)
 
subplot(3,2,4)
plot(ADEM.pU.x{1}(1,:),ADEM.pU.x{1}(2,:),ADEM.qU.x{1}(1,:),ADEM.qU.x{1}(2,:),':')
axis square
title('under tight priors','FontSize',16)
xlabel('displacement','FontSize',14)
axis([-1 1 -1 1]*2)
 
subplot(3,1,1)
hold on
plot(ADEM.pU.x{1}(1,:),ADEM.pU.x{1}(2,:),ADEM.qU.x{1}(1,:),ADEM.qU.x{1}(2,:),':')
hold off
 
subplot(3,2,5)
plot([1:N],DEM.qU.a{2},[1:N],DEM.qU.v{2},':',[1:N],DEM.pU.v{2},'-.')
axis square
title('action and cause','FontSize',16)
xlabel('time','FontSize',14)
legend({'action', 'perceived perturbation', 'true perturbation'})
axis([1 N -1 1])
 
subplot(3,2,6)
plot([1:N],ADEM.qU.a{2},[1:N],ADEM.qU.v{2},':',[1:N],ADEM.pU.v{2},'-.')
axis square
title('action and cause','FontSize',16)
xlabel('time','FontSize',14)
axis([1 N -1 1])

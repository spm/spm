% This demo re-visits the mountain car problem to show that adaptive
% (desired) behaviour can be prescribed in terms of loss-functions (i.e.
% reward functions of state-space).
% It exploits the fact that under the free-energy formulation, loss is
% divergence. This means that priors can be used to make certain parts of
% state-space costly (i.e. with high divergence) and others rewarding (low
% divergence). Active inference under these priors will lead to sampling of
% low cost states and (apparent) attractiveness of those states.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_mountaincar_loss.m 3757 2010-03-08 11:41:53Z guillaume $
 
% generative process (mountain car terrain)
%==========================================================================                        % switch for demo
clear
DEMO    = 0;
 
% range of position for later plotting
%--------------------------------------------------------------------------
dx      = 1/64;
x       = linspace(-2,2,1/dx);
xx      = x.^2;
 
% parameters of generative process
%--------------------------------------------------------------------------
P.a     = 0;
P.b     = [0 0];
P.c     = [0 0 0 0];
P.d     = 0;                                % action off
 
% level 1
%--------------------------------------------------------------------------
G(1).x  = [0; 0];
G(1).f  = inline('spm_fx_mountaincar(x,v,a,P)/2','x','v','a','P');
G(1).g  = inline('x','x','v','a','P');
G(1).pE = P;
G(1).V  = exp(16);                          % error precision
G(1).W  = exp(16);                          % error precision
 
% level 2
%--------------------------------------------------------------------------
G(2).a  = 0;                                % action
G(2).v  = 0;                                % inputs
G(2).V  = exp(16);
G       = spm_ADEM_M_set(G);
 

% generative model
%==========================================================================
clear P
 
% parameters and equations of motion
%--------------------------------------------------------------------------
P.x     =  0;
P.k     =  1;
P.p     =  -1/2;
P.q     =  0;
 
% level 1
%--------------------------------------------------------------------------
M(1).x  = [0; 0];
M(1).f  = inline('spm_mc_fx(x,v,P)/2','x','v','P');
M(1).g  = inline('x','x','v','P');
M(1).pE = P;
M(1).V  = exp(8);                           % error precision
M(1).W  = exp(6);                           % error precision
 
% level 2
%--------------------------------------------------------------------------
M(2).v  = 0;                                % inputs
M(2).V  = exp(16);
M       = spm_DEM_M_set(M);
 
 
% learn gradients with a flat loss-functions (priors on divergence)
%==========================================================================
DEM.U   = sparse(128,1);
DEM.C   = sparse(128,1);
DEM.G   = G;
DEM.M   = M;
 
if DEMO
    
    % enable learning by relaxing priors on parameters
    %----------------------------------------------------------------------
    DEM.M(1).pC = diag([ones(1,2) zeros(1,2)])/128;
    
    % initialise states (randomly) and integrate
    %----------------------------------------------------------------------
    for i = 1:8
        DEM.G(1).x  = 4*rand(2,1) - 2;
        DEM.M(1).x  = DEM.G(1).x;
        
        DEM         = spm_ADEM(DEM);
        DEM.M(1).pE = DEM.qP.P{1};
    end
else
    
    % use previously optimised parameters
    %----------------------------------------------------------------------
    DEM.M(1).pE.x = -0.2;
    DEM.M(1).pE.k =  0.6;
    
end
 
 
% plot results
%==========================================================================
DEM.G(1).x = [0; 1/2];
DEM.M(1).x = DEM.G(1).x;
DEM        = spm_ADEM(DEM);

spm_figure('GetWin','Figure 1'); 
spm_DEM_qU(DEM.qU)
 
% true and inferred position
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:)), hold on
plot(DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),':'),hold off
xlabel('position','Fontsize',14)
ylabel('velocity','Fontsize',14)
title('trajectories','Fontsize',16)
axis([-1 1 -1 1]*2)
axis square
 
% inferred potential
%--------------------------------------------------------------------------
dGdx = spm_mc_loss_G(x,DEM.qP.P{1});
 
% real potential
%--------------------------------------------------------------------------
dHdx = (x < 0).*(2*x + 1);
dHdx = (x > 0).*(1./(1 + 5*xx).^(1/2) - 5*xx./(1 + 5*xx).^(3/2) + (x/2).^4) + dHdx;
H    = cumsum(dHdx)*dx;
G    = cumsum(dGdx)*dx;
H    = H - min(H);
G    = G - min(G);
 
subplot(2,2,4)
plot(x,H,x,G,'-.')
xlabel('position','FontSize',14)
ylabel('height','FontSize',14)
title('implicit potential','FontSize',16)
axis square
drawnow
 
 
% enable action and cost-priors
%==========================================================================


% enable action (disable learning) and integrate
%--------------------------------------------------------------------------
DEM.G(1).pE.d = 1;
DEM.M(1).pE.p = 1/4;
DEM.M(1).pE.q = 32;
DEM.M(1).pC   = [];
DEM           = spm_ADEM(DEM);
 

% show trajectories
%==========================================================================

% loss-functions or priors
%--------------------------------------------------------------------------
c   = spm_mc_loss_C(x,DEM.M(1).pE);
C   = DEM.M(1).pE.p - DEM.M(1).pE.q*(1 - c);

% for plotting
%--------------------------------------------------------------------------
[Sx Sy Sz]   = sphere;
S{1}         = Sx/4 + 1;
S{2}         = Sy/4 + 0;
S{3}         = Sz/4 - 1;

% inference
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_DEM_qU(DEM.qU)
subplot(2,2,4),  hold on
plot(DEM.U,':'), hold off
 
% true and inferred position
%--------------------------------------------------------------------------
subplot(2,2,1)
% surf(S{1},S{2},S{3}), hold on
plot(1,0,'c.','MarkerSize',64), hold on
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:))
plot(DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),':'),hold off
xlabel('position','Fontsize',14)
ylabel('velocity','Fontsize',14)
title('trajectories','Fontsize',16)
axis([-1 1 -1 1]*2), shading interp
axis square
 
% cost-functions
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(x,C)
xlabel('position','FontSize',14)
ylabel('force','FontSize',14)
title('priors (cost-function)','FontSize',16)
axis square
drawnow
 

% and a few more examples
%--------------------------------------------------------------------------
for i = 1:4
 
    % active inference
    %----------------------------------------------------------------------
    DEM.G(1).x  = rand(2,1) - 1/2;
    DEM.M(1).x  = DEM.G(1).x;
    DEM         = spm_ADEM(DEM);
 
 
    % true and inferred position
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 2');
    subplot(2,2,1), hold on
    plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:),'Color',[0.8 0.8 1])

 
end
 
xlabel('position','Fontsize',14)
ylabel('velocity','Fontsize',14)
title('trajectories','Fontsize',16)
axis([-1 1 -1 1]*2)
axis square
 
% illustrate different behaviours under different precisions
%==========================================================================
 
% high and low exploration - attractiveness
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
DEM.G(1).x = [0; -1/2];
DEM.M(1).x = DEM.G(1).x;
 
W     = [2 4 6 8];
for i = 1:4
    
    % active inference
    %----------------------------------------------------------------------
    DEM.M(1).W = exp(W(i));
    DEM        = spm_ADEM(DEM);
    
    % true and inferred position
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 3');
    subplot(2,2,i)

    % surf(S{1},S{2},S{3}), hold on
    plot(1,0,'c.','MarkerSize',64), hold on
    plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:))
    plot(DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),':'),hold off
    xlabel('position','Fontsize',14)
    ylabel('velocity','Fontsize',14)
    title(sprintf('%s (%i)','trajectories',W(i)),'Fontsize',16)
    axis([-1 1 -1 1]*2), shading interp
    axis square

end


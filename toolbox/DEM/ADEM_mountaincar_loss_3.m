% This demo re-visits the mountain car problem to show that adaptive
% (desired) behaviour can be prescribed in terms of loss-functions (i.e.
% reward functions of state-space).
% It exploits the fact that under the free-energy formulation, loss is
% divergence. This means that priors can be used to make certain parts of
% state-space costly (i.e. with high divergence) and others rewarding (low
% divergence). Active inference under these priors will lead to sampling of
% low cost states and (apparent) attractiveness of those states.
%
% This is version three; where the cost function is inferred.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_mountaincar_loss_3.m 3333 2009-08-25 16:12:44Z karl $
 
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
np      = 8;
nq      = 16;
P.p     = sparse(1,np);
P.q     = sparse(1,nq);
fx      = inline('spm_mc_fx_3(x,v,P)','x','v','P');
gx      = inline('x.x','x','v','P');
x0.x    = [0; 0];
x0.c    = 0;
 
% level 1
%--------------------------------------------------------------------------
M(1).x  = x0;
M(1).f  = fx;
M(1).g  = gx;
M(1).pE = P;
M(1).V  = exp(8);                           % error precision
M(1).W  = diag(exp([8 4 16]));              % error precision
 
% level 2
%--------------------------------------------------------------------------
M(2).v  = 0;                                % inputs
M(2).V  = exp(16);
M       = spm_DEM_M_set(M);
 
 
% learn gradients with a flat loss-functions (priors on divergence)
%==========================================================================
U       = sparse(128,M(1).m);
DEM.U   = U;
DEM.C   = U;
DEM.G   = G;
DEM.M   = M;
 
if DEMO
    
    % enable learning by relaxing priors on parameters
    %----------------------------------------------------------------------
    DEM.M(1).pC = diag([ones(1,np) zeros(1,nq)])/128;
    
    % initialise states (randomly) and integrate
    %----------------------------------------------------------------------
    for i = 1:32
        DEM.G(1).x    = 4*rand(2,1) - 2;
        DEM.M(1).x.x  = DEM.G(1).x;
        
        DEM           = spm_ADEM(DEM);
        DEM.M(1).pE   = DEM.qP.P{1};
        
        disp(DEM.qP.P{1})
    end
else
    
    % use previously optimised parameters
    %----------------------------------------------------------------------
    DEM.M(1).pE.p = [2.7 1.7 0.74 -0.51 -0.85 0.08 -0.23 -1.15];
end
 
 
% plot results
%==========================================================================
spm_figure('GetWin','DEM');
 
DEM.G(1).x   = [0; 1/2];
DEM.M(1).x.x = DEM.G(1).x;
DEM          = spm_ADEM(DEM);
 
spm_DEM_qU(DEM.qU)
 
% true and inferred position
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:)), hold on
plot(DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),':'),hold off
xlabel('position','Fontsize',14)
ylabel('velcitiy','Fontsize',14)
title('trajectories','Fontsize',16)
axis([-1 1 -1 1]*2)
axis square
 
% inferred potential
%--------------------------------------------------------------------------
dGdx = spm_DEM_basis(x,[],DEM.qP.P{1}.p);
dGdx = -dGdx;
 
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
 
% loss-functions or priors
%--------------------------------------------------------------------------
q0    = sparse(nq,1) - 1;
q     = sparse(nq*3/4,1,-32,nq,1) + 1;
C0    = spm_DEM_basis(x,[],q0);
C     = spm_DEM_basis(x,[],q);
[i j] = min(C);
T     = x(j);
 
 
% enable action (disable learning) and integrate
%--------------------------------------------------------------------------
DEM.G(1).pE.d = 1;
DEM.M(1).pE.q = q;
DEM.M(1).pC   = [];
DEM           = spm_ADEM(DEM);
 
 
% show trajectories
%==========================================================================
 
% inference
%--------------------------------------------------------------------------
spm_figure('GetWin','Graphics');
spm_DEM_qU(DEM.qU)
subplot(2,2,4),  hold on
plot(U,':'), hold off
 
% true and inferred position
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:),T,0,'r.','Markersize',32), hold on
plot(DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),':'),hold off
xlabel('position','Fontsize',14)
ylabel('velcitiy','Fontsize',14)
title('trajectories','Fontsize',16)
axis([-1 1 -1 1]*2)
axis square
 
% cost-functions
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(x,C,x,C0,'b--')
xlabel('position','FontSize',14)
ylabel('force','FontSize',14)
title('priors (loss-functions)','FontSize',16)
axis square
drawnow
 
return
 
 
% and a few more examples
%--------------------------------------------------------------------------
for i = 1:4
 
    % active inference
    %----------------------------------------------------------------------
    DEM.G(1).x    = 2*rand(2,1) - 1;
    DEM.M(1).x.x  = DEM.G(1).x;
    DEM           = spm_ADEM(DEM);
 
 
    % true and inferred position
    %----------------------------------------------------------------------
    spm_figure('GetWin','Graphics');
    subplot(2,2,1), hold on
    plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:),'Color',[0.8 0.8 1])
    xlabel('position','Fontsize',14)
    ylabel('velcitiy','Fontsize',14)
    title('trajectories','Fontsize',16)
    axis([-1 1 -1 1]*2)
    axis square
 
end
 
 
% illustrate different behaviours under different precisions
%==========================================================================
 
% high and low exploration - attractiveness
%--------------------------------------------------------------------------
spm_figure('GetWin','Graphics'); clf
DEM.G(1).x = [0; 1/2];
DEM.M(1).x.x = DEM.G(1).x;
 
W     = [1 2 3 4];
for i = 1:4
    
    % active inference
    %----------------------------------------------------------------------
    DEM.M(1).W = exp(W(i));
    DEM        = spm_ADEM(DEM);
 
    % true and inferred position
    %----------------------------------------------------------------------
    spm_figure('GetWin','Graphics');
    subplot(2,2,i)
    plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:),T,0,'r.','Markersize',32), hold on
    plot(DEM.qU.x{1}(1,:),DEM.qU.x{1}(2,:),':'),hold off
    xlabel('position','Fontsize',14)
    ylabel('velcitiy','Fontsize',14)
    title('trajectories','Fontsize',16)
    axis([-1 1 -1 1]*2)
    axis square
 
end

% This demo re-visits the mountain car problem to show that adaptive
% (desired) behaviour can be prescribed in terms of loss-functions (i.e.
% reward functions of state-space).
% It exploits the fact that under the free-energy formulation, loss is
% divergence. This means that priors can be used to make certain parts of
% state-space costly (i.e. with high divergence) and others rewarding (low
% divergence). Active inference under these priors will lead to sampling of
% low cost states and (apparent) attractiveness of those states.
%
% This is version four; that includes a drive state.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_mountaincar_loss_4.m 4626 2012-01-24 20:55:59Z karl $
 
% generative process (mountain car terrain)
%==========================================================================                        % switch for demo
 
% parameters of generative process
%--------------------------------------------------------------------------
P.a     = 0;
P.b     = [0 0];
P.c     = [0 0 0 0];
P.d     = 1;                                % action on
 
fx      = inline('spm_mc_fxa_4(x,v,a,P)','x','v','a','P');
gx      = inline('[x.x; x.v; x.d]','x','v','a','P');
x0.x    = 0;
x0.v    = 0;
x0.p    = 0;
x0.d    = 0;


% level 1
%--------------------------------------------------------------------------
G(1).x  = x0;
G(1).f  = fx;
G(1).g  = gx;
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
clear P x0
 
% parameters (previously learned) and equations of motion
%--------------------------------------------------------------------------
P       = [2.7 1.7 0.74 -0.51 -0.85 0.08 -0.23 -1.15];
np      = length(P);
fx      = inline('spm_mc_fx_4(x,v,P)','x','v','P');
gx      = inline('[x.x; x.v; x.d]','x','v','P');
x0.x    = 0;
x0.v    = 0;
x0.c    = 0;
x0.d    = 0;
 
% level 1
%--------------------------------------------------------------------------
M(1).x  = x0;
M(1).f  = fx;
M(1).g  = gx;
M(1).pE = P;
M(1).V  = exp(8);                           % error precision
M(1).W  = diag(exp([8 4 16 16]));           % error precision
 
% level 2
%--------------------------------------------------------------------------
M(2).v  = 0;                                % inputs
M(2).V  = exp(16);
M       = spm_DEM_M_set(M);
 
 
% learn gradients with a flat loss-functions (priors on divergence)
%==========================================================================
N       = 1600;
U       = sparse(N,M(1).m);
DEM.U   = U;
DEM.C   = U;
DEM.G   = G;
DEM.M   = M;
DEM     = spm_ADEM(DEM);
 

% show dynamics
%==========================================================================
 
% inference
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU)
 
% true and inferred position
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:)),hold on
plot(   1,0,'r.','Markersize',32),      hold on
plot(-1/2,0,'g.','Markersize',16),      hold off
xlabel('position','Fontsize',14)
ylabel('velcitiy','Fontsize',14)
title('trajectories','Fontsize',16)
axis([-1 1 -1 1]*3)
axis square

% true position
%--------------------------------------------------------------------------
subplot(2,2,3)
plot3(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:),1:N), hold on
plot3(   1,0,1:64:N,'r.','Markersize',8),     hold on
plot3(-1/2,0,1:64:N,'g.','Markersize',8),     hold off
xlabel('position','Fontsize',14)
ylabel('velocity','Fontsize',14)
zlabel('time','Fontsize',14)
title('trajectories','Fontsize',16)
axis([-2 2 -2 2 0 N])
axis square


% real states
%==========================================================================
spm_figure('GetWin','DEM');
spm_DEM_qU(DEM.pU)
 

subplot(2,2,1)
plot3(   1,0,1:1/8:2,'r.'),                                       hold on
plot3(-1/2,0,1:1/8:2,'g.'),                                       hold on
plot3(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:),DEM.pU.x{1}(4,:)),        hold off
xlabel('position','Fontsize',14)
ylabel('velocity','Fontsize',14)
zlabel('satiety','Fontsize',14)
title('trajectories','Fontsize',16)
axis([-2 2 -2 2 0 8])



% cost function (see spm_mc_fx_4.m)
%--------------------------------------------------------------------------
subplot(2,1,2)
x     = -2:1/64:2;
d     =  0:1/64:2;
for i = 1:length(x)
    for j = 1:length(d)
        D      = spm_phi((1 - d(j))*8);
        A      = 2 - 32*exp(-(x(i) - 1).^2*32);
        C(i,j) = A*D - 1;
    end
end

surf(d,x,C)
shading interp
xlabel('drive','Fontsize',14)
ylabel('position','Fontsize',14)
title('cost-function','Fontsize',16)
axis square

return

% NOTES for graphics
%--------------------------------------------------------------------------
clear C
r     =  0:1/64:1;
d     =  0:1/64:4;
for i = 1:length(r)
    for j = 1:length(d)
        D      = spm_phi((1 - d(j))*8);
        A      = 2 - 32*r(i);
        C(i,j) = A*D - 1;
    end
end

imagesc(d,r,C)
shading interp
xlabel('satiety','Fontsize',14)
ylabel('reward','Fontsize',14)
axis xy


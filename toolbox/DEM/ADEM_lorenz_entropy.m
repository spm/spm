function ADEM_lorenz_entropy
% This demo shows how structure can be instilled from an environment. It
% uses an agent that optimises its recognition density over successive
% epochs of exposure to an environment. This environment causes the agent
% to flow on a Lorenz attractor with random perturbations. As the agent
% learns the causal regularities in its environment, it is better able to
% predict them and act to oppose the random effect. The result is that it
% is more robust to random forces and therefore exhibits states with lower
% entropy. This routine takes several minutes to run.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_lorenz_entropy.m 4804 2012-07-26 13:14:18Z karl $
 
% generative process (environment)
%==========================================================================
 
% set up
%--------------------------------------------------------------------------
G(1).E.s = 1/4;                             % smoothness
G(1).E.n = 4;                               % embedding
G(1).E.d = 2;                               % embedding
 
 
% dynamics
%--------------------------------------------------------------------------
fG      = '[v + a; 0; 0] + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64';
fM      = '[v    ; 0; 0] + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64';
 
% parameters
%--------------------------------------------------------------------------
PG      = [10; -8/3; 32];                   % Target parameters
PM      = [10;    0;  0];                   % initial parameters
 
% level 1
%--------------------------------------------------------------------------
G(1).x  = [1; 1; 24];
G(1).f  = inline(fG ,'x','v','a','P');
G(1).g  = inline('x','x','v','a','P');
G(1).pE = PG;
G(1).V  = exp(8);                           % error precision
G(1).W  = exp(16);                          % error precision
 
% level 2
%--------------------------------------------------------------------------
G(2).a  = 0;                                % action variables
G(2).v  = 0;                                % inputs
G(2).V  = exp(16);
G       = spm_ADEM_M_set(G);
 
 
% % plot flow fields and equilibrium densities
% %==========================================================================
% spm_figure('GetWin','Figure 1');
%
% x{1}    = linspace(-20,20,32);
% x{2}    = linspace(-32,32,32);
% x{3}    = linspace(  10,40,8);
%
%
% % controlled flow (P0)
% %--------------------------------------------------------------------------
% subplot(3,2,1)
% spm_fp_display_density(G,x);
% xlabel('position','Fontsize',12)
% ylabel('velocity','Fontsize',12)
% title('controlled','Fontsize',16)
 
 
% recognition model: learn the controlled environmental dynamics
%==========================================================================
 
% make a naive model (M)
%--------------------------------------------------------------------------
M       = G;
M(1).f  = inline(fM ,'x','v','P');
M(1).g  = inline('x','x','v','P');
M(1).pE = PM;
M(1).pC = eye(3)*4;
M(1).V  = exp(8);
M(1).W  = exp(8);
M       = spm_DEM_M_set(M);
 
% teach naive model by exposing it to a controlled environment (G)
%==========================================================================
 
% length of realization
%--------------------------------------------------------------------------
n     = 256;
U     = sparse(n,1);
 
DEM.M = M;
DEM.G = G;
DEM.C = U;
DEM.U = U;
 
% optimise recognition model in epochs
%--------------------------------------------------------------------------
DEM.M(1).E.nE = 1;
H     = sparse(128,128);
for i = 1:32
        
        % random perturbations
        %------------------------------------------------------------------
        DEM.C       = spm_conv(randn(n,1)*8,8);
        
        % integrate and update priors
        %------------------------------------------------------------------
        DEM         = spm_ADEM(DEM);
        DEM.M(1).pE = DEM.qP.P{1};
        DEM.M(1).x  = DEM.qU.x{1}(:,end);
        DEM.M(2).a  = DEM.qU.a{2}(:,end);
        
        % save free-energy and conditional expectations
        %------------------------------------------------------------------
        qP(:,i) = DEM.qP.P{1};
        F(i)    = DEM.F;
        
        % display states at the beginning and end
        %------------------------------------------------------------------
        spm_figure('GetWin','Figure 1');
  
        if i < 8
            subplot(3,2,3) 
            x = DEM.pU.v{1}(2,:);
            y = DEM.pU.v{1}(3,:);
            plot(x,y); hold on   
            title('states before','FontSize',16)
            a = axis;
            
        elseif i > 24
            subplot(3,2,4)
            x = DEM.pU.v{1}(2,:);
            y = DEM.pU.v{1}(3,:);
            H = H + sparse(64 + fix(x)',64 + fix(y)',1,128,128);
            plot(x,y); hold on   
            title('and after','FontSize',16)
            axis(a)
 
        end 
    end
 
 
 
% graphics
%==========================================================================
spm_figure('GetWin','Figure 1');
 
 
% plot free-energy
%--------------------------------------------------------------------------
t   = 1:length(F);
subplot(3,1,1)
plot(t,-F)
xlabel('time (epochs)','FontSize',12)
title('Free-energy','FontSize',16)
axis square
 
% and underlying conditional expectations
%--------------------------------------------------------------------------
subplot(3,2,5)
plot(t,qP,t,kron(ones(1,length(t)),PG),'-.')
xlabel('time (epochs)','FontSize',12)
title('conditional parameters','FontSize',16)
axis square
 
subplot(3,2,6)
imagesc(H')
xlabel('first state','FontSize',12)
ylabel('second state','FontSize',12)
title('ensemble density','FontSize',16)
axis square xy

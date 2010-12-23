% This demo shows how structure can be instilled from an environment. It 
% uses an agent that optimises its recognition density over successive 
% epochs of exposure to an environment. This environment causes the agent 
% to flow on a Lorenz attractor with random perturbations. As the agent 
% learns the causal regularities in its environment, it is better able to 
% predict them and act to oppose the random effect. The result is that it 
% is more robust to random forces and therefore exhibits states with lower 
% entropy.
%
% THIS ROUTINE IS NOT RUN IN DEMO MODE AND WILL TAKE A CONSIDEABLE TIME.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_lorenz_entropy.m 4146 2010-12-23 21:01:39Z karl $
 
% generative process (environment)
%==========================================================================
clear
DEMO     = 1;
G(1).E.s = 1/4;                        % smoothness
G(1).E.n = 6;                          % smoothness
G(1).E.d = 2;                          % smoothness
 
 
% dynamics
%--------------------------------------------------------------------------
fG      = '[v + a; 0; 0] + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64';
fM      = '[v    ; 0; 0] + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64';
 
% parameters
%--------------------------------------------------------------------------
PG      = [10; -8/3; 32];
PM      = [0;     0;  0];
 
% P(1): Prandtl number
% P(2): 8/3
% P(3): Rayleigh number
 
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
% spm_figure('GetWin','Graphics');
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
M(1).pC = eye(3)*2;
M(1).V  = [];
M(1).W  = [];
M(1).Q  = {speye(3,3)};
M(1).R  = {speye(3,3)};
M(1).hE = 0;
M(1).gE = 0;
M       = spm_DEM_M_set(M);
 
% teach naive model by exposing it to a controlled environment (G)
%--------------------------------------------------------------------------
clear DEM
 
% length of realization
%--------------------------------------------------------------------------
n     = 256;
U     = sparse(n,1);
 
M(1).E.nE = 1;
M(1).E.nM = 8;
 
DEM{1}.M = M;
DEM{1}.G = G;
DEM{1}.C = U;
DEM{1}.U = U;
 
% optimise recognition model in epochs
%--------------------------------------------------------------------------
if DEMO
    
    for i = 1:128
        
        % random perturbations
        %------------------------------------------------------------------
        DEM{i}.C           = spm_conv(randn(n,1)*8,8);
        
        % integrate and update priors
        %------------------------------------------------------------------
        DEM{i + 1}         = spm_ADEM(DEM{i});
        DEM{i + 1}.M(1).pE = DEM{i + 1}.qP.P{1};
        DEM{i + 1}.M(1).hE = DEM{i + 1}.qH.h{1};
        DEM{i + 1}.M(1).gE = DEM{i + 1}.qH.g{1};
        DEM{i + 1}.M(1).x  = DEM{i + 1}.qU.x{1}(:,end);
        DEM{i + 1}.M(2).a  = DEM{i + 1}.qU.a{2}(:,end);
        
        % display
        %------------------------------------------------------------------
        spm_DEM_qU(DEM{i + 1}.qU,DEM{i + 1}.pU)
        
    end
    
    save DEM_lorenz_entropy DEM
else
    load DEM_lorenz_entropy
end
 
% remove initial DEM and extract free-energy and conditional expectations
%--------------------------------------------------------------------------
DEM   = DEM(2:end);
for i = 1:length(DEM)
    qG(i)   = DEM{i}.qH.g{1};
    qH(i)   = DEM{i}.qH.h{1};
    qP(:,i) = DEM{i}.qP.P{1};
    F(i)    = DEM{i}.F;
end
 
% graphics
%==========================================================================
spm_figure('GetWin','Graphics');
t   = 1:length(DEM);
 
% plot free-energy
%--------------------------------------------------------------------------
subplot(3,1,1)
plot(t,-F)
xlabel('time (epochs)','FontSize',12)
title('Free-energy','FontSize',16)
axis square

% states
%--------------------------------------------------------------------------
subplot(3,2,3)
j = 1;
n = 16;
H = sparse(128,128);
for i = 1:n
    x = DEM{j + i}.pU.v{1}(2,:);
    y = DEM{j + i}.pU.v{1}(3,:);
    H = H + sparse(64 + fix(x)',64 + fix(y)',1,128,128);
    plot(x,y); hold on
end
hold off
title('states before','FontSize',16)
a = axis;

subplot(3,2,4)
j = length(t) - n;
H = sparse(128,128);
for i = 1:n
    x = DEM{j + i}.pU.v{1}(2,:);
    y = DEM{j + i}.pU.v{1}(3,:);
    H = H + sparse(64 + fix(x)',64 + fix(y)',1,128,128);
    plot(x,y); hold on
end
hold off
title('and after','FontSize',16)
axis(a)

 
% and underlying conditional expectations
%--------------------------------------------------------------------------
subplot(3,2,5)
plot(t,qP,t,kron(ones(1,length(t)),PG),'-.')
xlabel('time (epochs)','FontSize',12)
title('conditional parameters','FontSize',16)
axis square
 
subplot(3,2,6)
plot(t,qH,'-.',t,qG)
xlabel('time (epochs)','FontSize',12)
title('hyperparameters','FontSize',16)
axis square

% This demo illustrates the use of Lotka-Volterra form SHCs (Stable
% heteroclinic channel) to prescribe active sampling (inference). In this
% example each (unstable) fixed point in the SHC attracts the agent to
% points on the circumference of a circle.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_SHC_demo.m 3901 2010-05-27 16:14:36Z karl $
 
% generative process
%==========================================================================                        % switch for demo
clear
spm_figure('GetWin','Graphics');
 
 
% level 1
%--------------------------------------------------------------------------
G(1).x  = [1/2;0];
G(1).f  = inline('a','x','v','a','P');
G(1).g  = inline('x','x','v','a','P');
G(1).V  = exp(8);                           % error precision
G(1).W  = exp(8);                           % error precision
 
% level 2
%--------------------------------------------------------------------------
G(2).v  = 0;                                % inputs
G(2).a  = [0;0];                            % action
G(2).V  = exp(16);
G       = spm_ADEM_M_set(G);
 
 
% generative model
%==========================================================================                        
fx      = inline('spm_lotka_volterra(x,v,P)','x','v','P');
gx      = inline('spm_gx_SHC(x,v,P)','x','v','P');
gx      = inline('x.v','x','v','P');
 
% positions associated with each state (on unit circle)
%--------------------------------------------------------------------------
nx        = 8;
P.g(1,:)  = cos(2*pi*([1:nx]' - 1)/nx);
P.g(2,:)  = sin(2*pi*([1:nx]' - 1)/nx);
 
% parameters of succession
%--------------------------------------------------------------------------
P.f       =  spm_speye(nx,nx,-1) - spm_speye(nx,nx,+1);
P.f(nx,1) = -1;
P.f(1,nx) = +1;
P.f       =  P.f/2 - 1 + speye(nx);
 
% level 1
%--------------------------------------------------------------------------
M(1).x.x  = sparse(1,1,8,nx,1) - 8;
M(1).x.v  = [1; 0];
M(1).f    = fx;
M(1).g    = gx;
M(1).pE   = P;
M(1).V    = exp(4);                           % error precision
M(1).W    = exp(8);                           % error precision
 
% level 2
%--------------------------------------------------------------------------
M(2).v    = 0;                                % inputs
M(2).V    = exp(16);
M         = spm_DEM_M_set(M);
 
 
% ADEM
%==========================================================================
U       = sparse(128,M(1).m);
DEM.U   = U;
DEM.C   = U;
DEM.G   = G;
DEM.M   = M;
 
DEM     = spm_ADEM(DEM);
 
subplot(2,2,3)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:))
axis([-1 1 -1 1])
axis square

function ADEM_SHC_demo
% This demo illustrates the use of Lotka-Volterra form SHCs (Stable
% heteroclinic channel) to prescribe active sampling (inference). In this
% example each (unstable) fixed point in the SHC attracts the agent to
% points on the circumference of a circle.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_SHC_demo.m 6290 2014-12-20 22:11:50Z karl $
 
% generative process
%==========================================================================
M(1).E.n = 2;
M(1).E.d = 2;
M(1).E.s = 1;
 
% level 1
%--------------------------------------------------------------------------
G(1).x  = [1/2;0];
G(1).f  = @(x,v,a,P) a;
G(1).g  = @(x,v,a,P) x;
G(1).V  = exp(8);                           % error precision
G(1).W  = exp(8);                           % error precision
G(1).U  = 1;                                % error precision


% level 2
%--------------------------------------------------------------------------
G(2).v  = 0;                                % inputs
G(2).a  = [0;0];                            % action
G(2).V  = exp(16);

% generative model
%==========================================================================                        
 
% positions associated with each state (on unit circle)
%--------------------------------------------------------------------------
nx       = 8;
P.g(1,:) = cos(2*pi*((1:nx)' - 1)/nx);
P.g(2,:) = sin(2*pi*((1:nx)' - 1)/nx);
 
% parameters of succession
%--------------------------------------------------------------------------
P.f      =  spm_lotka_volterra(nx,1/2);
 
% level 1
%--------------------------------------------------------------------------
M(1).x   = sparse(1,1,8,nx,1) - 8;
M(1).f   = @(x,v,P) spm_lotka_volterra(x,v,P.f);
M(1).g   = @(x,v,P) P.g*spm_softmax(x);
M(1).pE  = P;
M(1).V   = exp(2);                           % error precision
M(1).W   = exp(8);                           % error precision
 
% level 2
%--------------------------------------------------------------------------
M(2).v = 0;                                   % inputs
M(2).V = exp(16);

 
% ADEM
%==========================================================================
U      = sparse(128,1);
DEM.U  = U;
DEM.C  = U;
DEM.G  = G;
DEM.M  = M;
DEM    = spm_ADEM(DEM);
 

% trajectory
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(DEM.pU.x{1}(1,:),DEM.pU.x{1}(2,:))
axis([-1 1 -1 1]*3/2)
axis square

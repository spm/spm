% Demo for phase-space reduction (cf. Centre-Manifold theory): In this
% example we show that Generalised filtering can be used to approximate the
% dynamics on the Centre-Manifold of a chaotic system. This example uses a
% nominal polynomial for the reduced dynamics (cf, order parameters) and
% exploits parameter optimisation to learn their coefficients. The approach
% appeals to Centre Manifold Theory and the enslaving principle by assuming
% that stable modes collapse sufficiently quickly to be treated as random
% fluctuations. This means that the unstable (slow modes or patterns)
% dynamics are sufficient to characterise the emergent order (e.g. pattern
% formation) and can be accessed using a generative model that tries to
% explain observed motion in terms of a few hidden states (i.e. order
% parameters) that implicitly enslave random fluctuations.
%
% THIS DEMO IS UNDER CONSTRUCTION
 
 
% non-hierarchical non-linear generative model (dynamic & chaotic)
%==========================================================================
spm_figure('GetWin','Figure 1');

% get model parameters and initial states
%--------------------------------------------------------------------------
x      = [4 4 28]';
P.A    = [-10  10  0;
           32 -1   0;
           0   0  -8/3];
P.B{1} = [ 0  0  0;
           0  0 -1;
           0  1  0];
P.B{2} = [ 0  0  0;
           0  0  0;
           0  0  0];
P.B{3} = [ 0  0  0;
           0  0  0;
           0  0  0];
P.C    = {};
 
V.A    = [ 1  1  1;
           1  1  1;
           1  1  1];
V.B{1} = [ 1  1  1;
           1  1  1;
           1  1  1];
V.B{2} = [ 0  1  1;
           0  1  1;
           0  1  1];
V.B{3} = [ 0  0  1;
           0  0  1;
           0  0  1];
V.C    = {};
pE     = spm_unvec(spm_vec(P)/2,P);
pC     = diag(spm_vec(V));
 
% level 1
%--------------------------------------------------------------------------
M(1).f  = inline('spm_fx_poly(x,v,P)/32','x','v','P');
M(1).g  = inline('x','x','v','P');
M(1).x  = x;
M(1).pE = P;
M(1).V  = exp(8);
M(1).W  = exp(8);
 
% level 2
%--------------------------------------------------------------------------
M(2).v  = 0;
M(2).V  = exp(16);
 
% create data
%==========================================================================
 
% create innovations & add causes
%--------------------------------------------------------------------------
N       = 128;
U       = sparse(1,N);
DEM     = spm_DEM_generate(M,U);
spm_DEM_qU(DEM.pU)
 
 
% DEM estimation
%==========================================================================
DEM.M(1).E.n  = 4;
DEM.M(1).E.s  = 1/8;
DEM.M(1).pE   = pE;
DEM.M(1).pC   = pC*exp(0);
DEM           = spm_DEM(DEM);
 
M(1).pE = DEM.qP.P{1};
SIM     = spm_DEM_generate(M,U);
 
 
% graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); 
spm_DEM_qU(SIM.pU)
 
subplot(2,2,3)
x      = DEM.pU.x{1};
plot(x(1,:),x(2,:))
title('true manifold dynamics')
xlabel('time')
axis square
 
subplot(2,2,4)
x      = SIM.pU.x{1};
plot(x(1,:),x(2,:))
title('estimated dynamics')
xlabel('time')
axis square
 
 
return
 
% Notes for creation of a globally coupled map
%==========================================================================
 
 
% Demo of synchronization manifold using coupled Lorenz attractors
%--------------------------------------------------------------------------
W    = 2;                            % amplitude of random fluctuations
N    = 16;                           % number of (Lorenz) oscillators
T    = 512;                          % number of time bins
dt   = 1/32;                         % time interval
 
% parameters (set s.d. P.t to 1/4 to see oscillator death)
%--------------------------------------------------------------------------
P.t  = randn(N,1)/16;                 % variations in log-rate constants
P.k  = 2;                            % global coupling parameter
 
% states
%--------------------------------------------------------------------------
x      = randn(3,N)*8;               % microstates
x(1,:) = x(1,:) + 0;
x(2,:) = x(2,:) + 0;
x(3,:) = x(3,:) + 28;
v      = 0;
 
 
% integrate
%--------------------------------------------------------------------------
for i = 1:T
    
    [dfdx f]  = spm_diff('spm_lorenz_k',x,v,P,1);
    f         = f + randn(N*3,1)*W;
    dx        = spm_dx(dfdx,f,dt);
    x         = x + spm_unvec(dx,x);
    y(:,i)    = x(:);               % microstates
    X(:,i)    = mean(x,2);          % macrostates
 
end
 
% plot
%--------------------------------------------------------------------------
subplot(2,1,1)
plot(y',':'), hold on
plot(X','k'), hold off
axis tight
 
t = T/4:T;
subplot(2,1,2)
plot(X(1,t),X(2,t),'k')
axis square
axis([-16 16 -16 16])
 
return
 
% plot synchronisation manifold
%--------------------------------------------------------------------------
clf
subplot(2,1,1)
for i = 1:12
    plot(y(i*3 - 2,t),y(i*3 + 1,t),'.','Color',[1 1/2 1/2]),hold on
end
plot(X(1,t),X(1,t),'.k'),hold off
axis square
axis([-16 16 -16 16])
 
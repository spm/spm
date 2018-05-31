function DEM_coupled_oscillators
% Dual estimation of the Lorenz system: Cross-validation of Laplace schemes
%__________________________________________________________________________
% Inversion of the Lorenz attractor with DEM, LAP and SCKS schemes: This
% demo tackles the difficult problem of deconvolving (chaotic) hidden states
% from a single response variable, while estimating the parameters of the
% underlying equations of motion. It calls generalised filtering, DEM and
% a state-of-the-art Bayesian smoother (SCKS).  This example is chosen to
% show that it is, in principle, possible to perform dual estimation in the
% context of chaotic dynamics (although small variations in this problem
% will cause the schemes to fail due it its inherently nonlinear nature and
% non-identifiability); however, the results are imperfect.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_coupled_oscillators.m 7322 2018-05-31 09:47:15Z karl $
 
 
% specify states and parameters
%==========================================================================
n     = 2;                               % number of sources (oscillators)
x.r   = zeros(n,1);                      % amplitude
x.p   = zeros(n,1);                      % phase differential
x.w   = zeros(1,1);                      % phase

% model parameters
%--------------------------------------------------------------------------
P.L   = eye(n,n);                        % lead field (measurement mapping)
P.Ap  = ones(n,n)/64 - eye(n,n)/16;          % amplitude coupling
P.Ar  = spm_speye(n,n,-1); 
P.Ar  = P.Ar/16 + P.Ar'/4;        % phase coupling
P.C   = [1; 0]/32;                          % exogenous input
P.r   = 1/32;                             % weak amplitude
P.w   = 2*pi/32;                          % intrinsic frequency

% % observation function
%--------------------------------------------------------------------------
g = @(x,v,P) P.L*((x.r).*cos(x.p + x.w));
g = @(x,v,P) [x.r; x.p];

% equations of motion
%--------------------------------------------------------------------------
f = @(x,v,P) [P.Ap*(x.r - P.r) + P.C*v;
              sum(P.Ar.*sin(bsxfun(@minus,x.p,x.p')),2) - P.C*v; ...
              P.w];

% causes or exogenous input
%--------------------------------------------------------------------------
N = 256;                                 % number of time points
U = exp(-((1:N) - N/4).^2/(2*(N/8)^2));             % exogenous input


E.n     = 4;
E.d     = 1;
E.nN    = 8;
E.s     = 1;

% first level dynamics
%--------------------------------------------------------------------------
M(1).E  = E;
M(1).x  = x;
M(1).f  = f;
M(1).g  = g;
M(1).pE = P;
M(1).V  = exp(12);
M(1).W  = exp(16);

% second level – causes or exogenous forcing term
%--------------------------------------------------------------------------
M(2).v = 0;
M(2).V = exp(16);

% create data
%==========================================================================
DEM    = spm_DEM_generate(M,U,P);
 
% show data
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.pU);

T  = (1:N);

subplot(4,2,2), plot(T,DEM.pU.x{1}(1:2,:)')
title('hidden amplitude','FontSize',16)
xlabel('time (seconds)'), spm_axis tight, box off

subplot(4,2,4), plot(T,DEM.pU.x{1}(3:4,:)')
title('phase difference','FontSize',16)
xlabel('time (seconds)'), spm_axis tight, box off

subplot(4,2,6), plot(T,abs(hilbert(full(DEM.Y)')))
title('response amplitude','FontSize',16)
xlabel('time (seconds)'), spm_axis tight, box off

subplot(4,2,8), plot(T,unwrap(angle(hilbert(full(DEM.Y)'))))
title('unwrapped phase','FontSize',16)
xlabel('time (seconds)'), spm_axis tight, box off
drawnow

% initialization of parameters (True values: pE =[18;-4;46.92])
%--------------------------------------------------------------------------
pE       = P;                          % prior parameters
pC       = spm_zeros(P);               % prior variance 

pE.Ar    = zeros(n,n);
pE.Ap    = -speye(n,n)/16;
pC.Ar    = (P.Ar ~= 0);
pC.Ap    = (P.Ap ~= 0);

DEM.M(1).pE = pE;
DEM.M(1).pC = sparse(diag(spm_vec(pC)));
DEM.M(1).V  = exp(4);
DEM.M(1).W  = diag(exp([4 4 4 4 32]));

DEM.U  = U;
 
 
% Inversion: 
%==========================================================================
LAP     = spm_LAP(DEM);

% Show estimates of states
%--------------------------------------------------------------------------
spm_figure('GetWin','spm_LAP'); spm_DEM_qU(LAP.qU,LAP.pU)


% and parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Parameters'); spm_DEM_qP(LAP.qP,LAP.pP)
subplot(2,1,1),legend('mean','90% CI')
title('Estimated and true (black) parameters','FontSize',16)



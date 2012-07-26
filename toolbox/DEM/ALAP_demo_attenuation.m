function DEM_demo_ALAP
% This demonstration is essentially the same as DEM_demo_LAP – however
% here, we compare two generalised filtering schemes that are implemented
% very differently: the first integrates the generative process in
% parallel with the inversion, while the standard spm_LAP scheme inverts a
% model given pre-generated data. The advantage of generating and modelling
% data  contemporaneously is that it allows the inversion scheme to couple
% back to the generative process through action (see active inference
% schemes): spm_ALAP.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ALAP_demo_attenuation.m 4804 2012-07-26 13:14:18Z karl $
 
% process (G) and model (M)
%==========================================================================

% set dimensions for generalised coordinates
%--------------------------------------------------------------------------
G(1).E.d        = 2;                   % approximation order
G(1).E.n        = 4;                   % embedding order
G(1).E.s        = 1/2;                 % temporal smoothness

M(1).E.d        = 2;                   % approximation order
M(1).E.n        = 4;                   % embedding order
M(1).E.s        = 1/2;                 % temporal smoothness
M(1).E.method.x = 1;                   % state-dependent noise
 
G(1).f  = inline('[a; v] - x/4','x','v','a','P');
G(1).g  = inline('[x(1); sum(x)]','x','v','a','P');
G(1).x  = [0; 0];                      % hidden state
G(1).v  = [0; 0];                      % hidden cause
G(1).V  = exp(8);                      % precision (noise)
G(1).W  = exp(8);                      % precision (states)
 
 
% level 2; causes
%--------------------------------------------------------------------------
G(2).v  = 0;                           % hidden cause
G(2).V  = exp(16);
 
 
% state-dependent precision (attentional bias) in generative model (M):
%--------------------------------------------------------------------------
M(1).f  = inline('v - x/4','x','v','P');
M(1).g  = inline('[x(1); sum(x)]','x','v','P');
M(1).x  = [0; 0];                      % hidden state
M(1).v  = [0; 0];                      % hidden cause
M(1).W  = exp(8);                      % precision (states)
M(1).ph = inline('[1; 1]*(8 + h*x(1))','x','v','h','M');


% level 2; causes
%--------------------------------------------------------------------------
M(2).v  = [0; 0];                      % hidden cause
M(2).V  = exp(4);


% free hyperparameters
%--------------------------------------------------------------------------
M(1).hE = 8;

 
% hidden cause and prior expectations
%========================================================================== 
N      = 32;
U      = exp(-((1:N) - 12).^2/(2.^2));

% invert
%==========================================================================
DEM.M  = M;
DEM.G  = G;
DEM.C  = C;
DEM.U  = U;

% generate and filter responses
%-------------------------------------------------------------------------- 
LAP    = spm_ALAP(DEM);




 
% Show results for LAP (standard scheme)
%==========================================================================
spm_figure('GetWin','Figure 1: Generalised filtering - standard scheme');
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
 
% parameters
%--------------------------------------------------------------------------
qP    = spm_vec(DEM.qP.P);
qP    = qP(ip);
tP    = spm_vec(DEM.pP.P);
tP    = tP(ip);
 
subplot(2,2,4)
bar([tP qP])
axis square
legend('true','GF – standard')
title('parameters','FontSize',16)
 
cq    = 1.64*sqrt(diag(DEM.qP.C(ip,ip)));
for i = 1:length(qP),hold on
    plot([i i] + 1/8,qP(i) + [-1 1]*cq(i),'LineWidth',4,'color','r')
end, hold off
 
 
% Show results for ALAP (parallel scheme)
%==========================================================================
spm_figure('GetWin','Figure 2: Generalised filtering – parallel scheme');
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(LAP.qU,LAP.pU)
 
% parameters
%--------------------------------------------------------------------------
qP    = spm_vec(LAP.qP.P);
qP    = qP(ip);
tP    = spm_vec(LAP.pP.P);
tP    = tP(ip);
 
subplot(2,2,4)
bar([tP qP])
axis square
legend('true','GF – parallel')
title('parameters','FontSize',16)
 
cq    = 1.64*sqrt(diag(LAP.qP.C(ip,ip)));
for i = 1:length(qP),hold on
    plot([i i] + 1/8,qP(i) + [-1 1]*cq(i),'LineWidth',4,'color','r')
end, hold off
 
% Compare
%==========================================================================
spm_figure('GetWin','Figure 3: Comparison of integration schemes');
 
% hyperparameters
%--------------------------------------------------------------------------
qL    = spm_vec({LAP.qH.h LAP.qH.g});
qD    = spm_vec({DEM.qH.h DEM.qH.g});
vL    = spm_vec({LAP.qH.V LAP.qH.W});
vD    = spm_vec({DEM.qH.V DEM.qH.W});
qh    = log([G(1).V; G(1).W]);
 
 
subplot(2,2,1)
bar([qh qL qD])
axis square
legend('true','parallel','standard')
title('log-precisions','FontSize',16)
 
cq    = 1.64*sqrt(vL);
for i = 1:length(qL),hold on
    plot([i i] + 0,qL(i) + [-1 1]*cq(i),'LineWidth',4,'color','r')
end, hold off
 
cq    = 1.64*sqrt(vD);
for i = 1:length(qD),hold on
    plot([i i] + 1/4,qD(i) + [-1 1]*cq(i),'LineWidth',4,'color','r')
end, hold off
 
% Log-evidence
%--------------------------------------------------------------------------
subplot(2,2,2)
nL   = length(LAP.F);
nD   = length(DEM.F);
plot(1:nL,LAP.F,1:nD,DEM.F)
axis square
legend('parallel (F)','standard (F)')
title('log-evidence ','FontSize',16)
xlabel('iteration','FontSize',12)

function DEM_demo_ALAP
% This demonstration is essentially the same as DEM_demo_LAP - however
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
% $Id: DEM_demo_ALAP.m 7679 2019-10-24 15:54:07Z spm $
 
% get basic convolution model
%==========================================================================
M       = spm_DEM_M('convolution model');

% gradient functions for speed (not implemented here)
%--------------------------------------------------------------------------
% M(1).fx = inline('P.f','x','v','P');
% M(1).fv = inline('P.h','x','v','P');
% M(1).gx = inline('P.g','x','v','P');
% M(1).gv = inline('sparse(4,1)','x','v','P');

M(1).E.nN = 8;                                 % number of E steps
M(1).E.nE = 8;                                 % number of E steps
M(1).E.nD = 1;                                 % number of time steps
M(1).E.s  = 1;                                 % smoothness
G(1).E.s  = 1;                                 % smoothness
M(1).E.d  = 2;                                 % order
M(1).E.n  = 6;                                 % order

 
% free parameters
%--------------------------------------------------------------------------
P       = M(1).pE;                             % true parameters
ip      = [1 2 5 9];                           % free parameters
pE      = spm_vec(P);
np      = length(pE);
pE(ip)  = 0;
pE      = spm_unvec(pE,P);
pC      = sparse(ip,ip,exp(4),np,np);
M(1).pE = pE;
M(1).pC = pC;
 
% free hyperparameters
%--------------------------------------------------------------------------
M(1).Q  = {speye(M(1).l,M(1).l)};
M(1).R  = {speye(M(1).n,M(1).n)};
M(1).hE = 8;
M(1).gE = 6;
M(1).hC = 1/4;
M(1).gC = 1/4;
 
% generative process
%==========================================================================
G(1).f  = M(1).f;
G(1).g  = M(1).g;
G(1).x  = M(1).x;
G(1).V  = exp(8);
G(1).W  = exp(6);
G(1).pE = P;

G(2).v  = 0;
G(2).V  = exp(16);


% hidden cause
%-------------------------------------------------------------------------- 
N      = 32;
U      = exp(-((1:N) - 12).^2/(2.^2));

% invert
%==========================================================================
DEM.M  = M;
DEM.G  = G;
DEM.C  = U;

% generate and filter responses
%-------------------------------------------------------------------------- 
LAP    = spm_ALAP(DEM);

% filter generated responses
%-------------------------------------------------------------------------- 
DEM.Y  = LAP.Y;
DEM.pU = LAP.pU;
DEM.pP = LAP.pP;
DEM    = spm_LAP(DEM);

 
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
legend('true','GF - standard')
title('parameters','FontSize',16)
 
cq    = 1.64*sqrt(diag(DEM.qP.C(ip,ip)));
for i = 1:length(qP),hold on
    plot([i i] + 1/8,qP(i) + [-1 1]*cq(i),'LineWidth',4,'color','r')
end, hold off
 
 
% Show results for ALAP (parallel scheme)
%==========================================================================
spm_figure('GetWin','Figure 2: Generalised filtering - parallel scheme');
 
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
legend('true','GF - parallel')
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

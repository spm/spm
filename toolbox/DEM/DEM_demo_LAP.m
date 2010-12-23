% This demonstration compares Generalised filtering under the Laplace
% assumption (spm_LAP) with variational filtering under the same fixed form
% approximation (i.e. DEM). We use a simple linear convolution model to
% illustrate the differences and similarities. The key difference between
% the two schemes lies (in this example) lies in estimates of conditional
% uncertainty. spm_LAP is must less over-confident because it eschews the
% means-field approximation implicit in DEM. The demonstration addresses 
% quadruple estimation of hidden states, exogenous input, parameters and 
% log-precisions (and, for spm_LAP, log-smoothness)
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_LAP.m 4146 2010-12-23 21:01:39Z karl $
 
% get basic convolution model
%==========================================================================
M       = spm_DEM_M('convolution model');
 
% gradient functions for speed (not implemented here)
%--------------------------------------------------------------------------
% M(1).fx = inline('P.f','x','v','P');
% M(1).fv = inline('P.h','x','v','P');
% M(1).gx = inline('P.g','x','v','P');
% M(1).gv = inline('sparse(4,1)','x','v','P');
 
% free parameters
%--------------------------------------------------------------------------
P       = M(1).pE;                            % true parameters
ip      = [1 2 5 9];                          % free parameters
pE      = spm_vec(P);
np      = length(pE);
pE(ip)  = 0;
pE      = spm_unvec(pE,P);
pC      = sparse(ip,ip,exp(8),np,np);
M(1).pE = pE;
M(1).pC = pC;
 
% free hyperparameters+
%--------------------------------------------------------------------------
M(1).Q  = {speye(M(1).l,M(1).l)};
M(1).R  = {speye(M(1).n,M(1).n)};
M(1).hE = 6;
M(1).gE = 6;
M(1).hC = 1/4;
M(1).gC = 1/4;
 
% generate data and invert
%==========================================================================
M(1).E.nN = 16;                                % number of E steps
M(1).E.nE = 16;                                % number of E steps
M(1).E.nD = 1;                                 % number of time steps
M(1).E.s  = 1;                                 % smoothness
M(1).E.d  = 2;                                 % order
M(1).E.n  = 6;                                 % order
 
N         = 32;                                % length of data sequence
U         = exp(-([1:N] - 12).^2/(2.^2));      % this is the Gaussian cause
DEM       = spm_DEM_generate(M,U,{P},{8,16},{6});
 
 
% invert
%==========================================================================
spm_figure('GetWin','DEM');
LAP       = spm_LAP(DEM);
DEM       = spm_DEM(DEM);
 
% Show results for DEM
%==========================================================================
spm_figure('GetWin','Figure 1: DEM - mean-field assumption');
 
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
legend('true','DEM')
title('parameters','FontSize',16)
 
cq    = 1.64*sqrt(diag(DEM.qP.C(ip,ip)));
for i = 1:length(qP),hold on
    plot([i i] + 1/8,qP(i) + [-1 1]*cq(i),'LineWidth',4,'color','r')
end, hold off
 
 
% Show results for LAP
%==========================================================================
spm_figure('GetWin','Figure 2: Generalised filtering (GF)');
 
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
legend('true','GF')
title('parameters','FontSize',16)
 
cq    = 1.64*sqrt(diag(LAP.qP.C(ip,ip)));
for i = 1:length(qP),hold on
    plot([i i] + 1/8,qP(i) + [-1 1]*cq(i),'LineWidth',4,'color','r')
end, hold off
 
% Compare
%==========================================================================
spm_figure('GetWin','Figure 3: Comparison of DEM and GF');
 
% hyperparameters
%--------------------------------------------------------------------------
qL    = spm_vec({LAP.qH.h LAP.qH.g});
qD    = spm_vec({DEM.qH.h DEM.qH.g});
vL    = spm_vec({LAP.qH.V LAP.qH.W});
vD    = spm_vec({DEM.qH.V DEM.qH.W});
qh    = spm_vec({DEM.pH.h{1} DEM.pH.g{1}});
 
 
subplot(2,2,1)
bar([qh qL qD])
axis square
legend('true','GF','DEM')
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
plot(1:nL,LAP.F,1:nD,DEM.F,1:nL,LAP.S,1:nD,DEM.S)
axis square
legend('GF (F)','DEM (F)','GF (S)','DEM(S)')
title('log-evidence ','FontSize',16)
xlabel('iteration','FontSize',12)

% Triple estimation of states, parameters and hyperparameters:
% This demo focuses estimating both the states and parameters to furnish a
% complete system identification, given only the form of the system and its
% responses to unknown input (c.f., DEM_demo_EM, which uses known inputs)
 
% get basic convolution model
%==========================================================================
M       = spm_DEM_M('convolution model');
 
% gradient functions for speed
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
pC      = sparse(ip,ip,32,np,np);
M(1).pE = pE;
M(1).pC = pC;
 
% free hyperparameters
%--------------------------------------------------------------------------
M(1).Q  = {speye(M(1).l,M(1).l)};
M(1).R  = {speye(M(1).n,M(1).n)};
M(1).hE = 2;
M(1).gE = 2;
M(1).hC = 8;
M(1).gC = 8;
M(1).V  = 0;
M(1).W  = 0;

% generate data and invert
%==========================================================================
M(1).E.nN = 16;                                % number of time steps
M(1).E.nD = 1;                                 % number of time steps
M(1).E.s  = 1/2;                               % smoothness
M(1).E.d  = 2;                                 % order
M(1).E.n  = 4;                                 % order
 
N         = 32;                                % length of data sequence
U         = exp(-([1:N] - 12).^2/(2.^2));      % this is the Gaussian cause
DEM       = spm_DEM_generate(M,U,{P},{8,16},{8});


% invert
%==========================================================================
LAP       = spm_LAP(DEM);
DEM       = spm_DEM(DEM);
 
% Show results for DEM
%==========================================================================
 
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
spm_figure('GetWin','SI');
 
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
legend('true','LAP')
title('parameters','FontSize',16)
 
cq    = 1.64*sqrt(diag(LAP.qP.C(ip,ip)));
for i = 1:length(qP),hold on
    plot([i i] + 1/8,qP(i) + [-1 1]*cq(i),'LineWidth',4,'color','r')
end, hold off
 
% Compare
%==========================================================================
spm_figure('GetWin','Laplace');
 
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
legend('true','LAP','DEM')
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
plot(1:length(LAP.F),LAP.F,1:length(DEM.F),DEM.F,'-.')
axis square
legend('LAP','DEM')
title('log-evidence ','FontSize',16)

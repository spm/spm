function DEM_demo_EM
% Dual estimation of parameters and hyperparameters; under known causes:
% This demo focuses on conditional parameter estimation with DEM and
% provides a comparative evaluation using EM.  This proceeds by removing
% uncertainly about the input so that the D-step can be discounted.

 
% get basic convolution model
%==========================================================================
M       = spm_DEM_M('convolution model');
 
% free parameters
%--------------------------------------------------------------------------
P       = M(1).pE;                            % true parameters
ip      = [2 5];                              % free parameters
pE      = spm_vec(P);
pE(ip)  = 0;
np      = length(pE);
pE      = spm_unvec(pE,P);
pC      = sparse(ip,ip,exp(8),np,np);
M(1).pE = pE;
M(1).pC = pC;
 
% free hyperparameters
%--------------------------------------------------------------------------
M(1).Q  = {speye(M(1).l,M(1).l)};
M(1).R  = {speye(M(1).n,M(1).n)};

% level 2
%--------------------------------------------------------------------------
M(2).l  = 1;                                  % inputs
M(2).V  = exp(16);                            % very precise causes
 

% and generate data
%==========================================================================
N       = 32;                                 % length of data sequence
U       = exp(-([1:N] - 12).^2/(2.^2));       % this is the Gaussian cause
DEM     = spm_DEM_generate(M,U,{P},{8,32},{32});


% invert model
%==========================================================================
DEM.U   = U;
DEM     = spm_DEM(DEM);

% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)


% EM: spm_nlsi_GN
%==========================================================================
G.f   =  inline('P.f*x + P.h*u','x','u','P','M');
G.g   =  inline('P.g*x','x','u','P','M');
G.m   =  DEM.M(1).m;
G.n   =  DEM.M(1).n;
G.l   =  DEM.M(1).l;
G.x   =  DEM.M(1).x;
G.pE  =  DEM.M(1).pE;
G.pC  =  DEM.M(1).pC;
G.hE  = -DEM.M(1).hE;
 
% exogenous inputs
%--------------------------------------------------------------------------
GU.u  = U';
GU.dt = 1;
 
% data and serial correlations
%--------------------------------------------------------------------------
t     = ((1:N) - 1);
K     = toeplitz(exp(-t.^2/(2*M(1).E.s^2)));
Q     = K*K';
 
GY.y  = DEM.Y';
GY.X0 = DEM.X';
GY.dt = 1;
GY.Q  = {kron(speye(G.l,G.l),Q)};
 
 
% EM with a Gauss-Newton-like optimization of free energy
%==========================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN(G,GU,GY);
 
% parameters
%--------------------------------------------------------------------------
ip    = [2 5];
qP    = spm_vec(DEM.qP.P);
qP    = qP(ip);
tP    = spm_vec(DEM.pP.P);
tP    = tP(ip);
pP    = spm_vec(DEM.M(1).pE);
pP    = pP(ip);
eP    = spm_vec(Ep);
eP    = eP(ip);
 
spm_figure('GetWin','DEM');
subplot(2,2,4)
bar([tP qP eP])
axis square
legend('true','DEM','EM')
title('parameters','FontSize',16)
 
cq    = 1.64*sqrt(diag(DEM.qP.C(ip,ip)));
ce    = 1.64*sqrt(diag(Cp(ip,ip)));
hold on
for i = 1:length(qP)
    plot([i i],       qP(i) + [-1 1]*cq(i),'LineWidth',8,'color','r')
    plot([i i] + 1/4, eP(i) + [-1 1]*ce(i),'LineWidth',8,'color','r')
end
hold off
 

return


% repeat for several realizations
%==========================================================================
clear QP EP QH EH
for i = 1:8
 
    % generate new data and DEM
    %----------------------------------------------------------------------
    DEM     = spm_DEM_generate(M,U,{P},{8,32},{32});
    DEM.U   = U;
    DEM     = spm_DEM(DEM);
 
    % EM
    %----------------------------------------------------------------------
    GY.y  = DEM.Y';
    [Ep,Cp,Eh,F] = spm_nlsi_GN(G,GU,GY);
 
    % retain parameter estimates
    %----------------------------------------------------------------------
    qP      = spm_vec(DEM.qP.P);
    qP      = qP(ip);
    eP      = spm_vec(Ep);
    eP      = eP(ip);
 
    QP(:,i) = qP;
    EP(:,i) = eP;
 
    QH(i) = DEM.qH.h{1}(1);
    EH(i) = Eh(1);
end
 
spm_figure('GetWin','Figure 1');

subplot(2,1,1)
bar(tP,'FaceColor',[1 1 1]*.9,'EdgeColor',[1 1 1]*.9)
hold on
plot([1 2] - 1/8,EP,'r.',[1 2] + 1/4,QP,'k.','Markersize',16)
hold off
axis square
set(gca,'XLim',[0 3])
legend('true','EM','DEM')
title('conditional estimates','FontSize',16)



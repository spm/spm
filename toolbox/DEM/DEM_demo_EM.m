% This demo focuses on conditional parameter estimation with DEM and
% provides a comparative evaluation using EM.  This proceeds by removing
% uncertainly about the input so that the D-step can be discounted.
%==========================================================================
clear M
 
% basic convolution model
%==========================================================================
 
% level 1
%--------------------------------------------------------------------------
M(1).E.linear = 1;

M(1).E.n  = 8;                                % embedding
M(1).E.d  = 3;                                % restriction
M(1).E.s  = 1/2;                              % smoothness
M(1).E.nE = 1;                                % E-steps
 
% level 1
%--------------------------------------------------------------------------
P.f     = [-1  4   ;                          % the Jacobian for the
           -2 -1]/4;                          % hidden sates
P.g     = [spm_dctmtx(4,2)]/4;                % the mixing parameters
P.h     = [1; 0];                             % input parameter
np      = length(spm_vec(P));
ng      = size(P.g,1);
nx      = size(P.f,1);
nv      = size(P.h,2);
 
% free parameters
%--------------------------------------------------------------------------
ip      = [2 5];                              % free parameters
pE      = spm_vec(P);
pE(ip)  = 0;
pE      = spm_unvec(pE,P);
pC      = sparse(ip,ip,exp(8),np,np);
 
M(1).pE = pE;
M(1).pC = pC;                                 % The prior covariance
 
M(1).n  = nx;
M(1).f  = inline('P.f*x + P.h*v','x','v','P');
M(1).g  = inline('P.g*x','x','v','P');
M(1).Q  = {speye(ng,ng)};
M(1).R  = {speye(nx,nx)};
 
 
% level 2
%--------------------------------------------------------------------------
M(2).l  = 1;                                  % inputs
M(2).V  = exp(16);
 
 
% and generate data
%==========================================================================
N       = 32;                                 % length of data sequence
U       = exp(-([1:N] - 12).^2/(2.^2));       % this is the Gaussian cause
DEM     = spm_DEM_generate(M,U,{P},{8,32},{32});
 
% display
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.pU)
 
 
% invert model
%==========================================================================
DEM.U   = U;
DEM     = spm_DEM(DEM);
 
 
% overlay true values
%--------------------------------------------------------------------------
subplot(2,2,2)
hold on
plot([1:N],DEM.pU.x{1},'linewidth',2,'color',[1 1 1]/2)
hold off
 
subplot(2,2,3)
hold on
plot([1:N],DEM.pU.v{2},'linewidth',2,'color',[1 1 1]/2)
hold off
 
drawnow
 
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
t     = ([1:N] - 1);
K     = toeplitz(exp(-t.^2/(2*M(1).E.s^2)));
Q     = K*K';
 
GY.y  = DEM.Y';
GY.X0 = DEM.X';
GY.dt = 1;
GY.Q  = {kron(speye(ng,ng),Q)};
 
 
% EM with a Gauss-Newton-like optimization of free energy
%==========================================================================
[Ep,Cp,S,F] = spm_nlsi_GN(G,GU,GY);
 
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
 
f = spm_figure;
subplot(2,1,1)
bar([tP qP eP])
axis square
legend('true','DEM','EM')
title('parameters')
 
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
    [Ep,Cp,S,F] = spm_nlsi_GN(G,GU,GY);
 
    % retain parameter estimates
    %----------------------------------------------------------------------
    qP      = spm_vec(DEM.qP.P);
    qP      = qP(ip);
    eP      = spm_vec(Ep);
    eP      = eP(ip);
 
    QP(:,i) = qP;
    EP(:,i) = eP;
 
    QH(i) = DEM.qH.h{1}(1);
    EH(i) = S(1);
end
 
figure(f)
subplot(2,1,2)
bar(tP,'FaceColor',[1 1 1]*.9,'EdgeColor',[1 1 1]*.9)
hold on
plot([1 2] - 1/8,EP,'r.',[1 2] + 1/4,QP,'k.','Markersize',16)
hold off
axis square
set(gca,'XLim',[0 3])
legend('true','EM','DEM')
title('conditional estimates')
 
%Triple estimation
%==========================================================================

% now repeat but without using the inputs
%--------------------------------------------------------------------------
M(2).V    = exp(-8);
M(1).E.nE = 1;
M(1).E.nN = 32;
DEM       = spm_DEM_generate(M,U,{P},{8,32},{32});
DEM       = spm_DEM(DEM);

% overlay true values
%--------------------------------------------------------------------------
subplot(2,2,2)
hold on
plot([1:N],DEM.pU.x{1},'linewidth',2,'color',[1 1 1]/2)
hold off
 
subplot(2,2,3)
hold on
plot([1:N],DEM.pU.v{2},'linewidth',2,'color',[1 1 1]/2)
hold off

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
title('parameters')
 
cq    = 1.64*sqrt(diag(DEM.qP.C(ip,ip)));
hold on
for i = 1:length(qP)
    plot([i i] + 1/8,qP(i) + [-1 1]*cq(i),'LineWidth',8,'color','r')
end
hold off


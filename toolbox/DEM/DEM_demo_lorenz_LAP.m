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
% $Id: DEM_demo_lorenz_LAP.m 4146 2010-12-23 21:01:39Z karl $
 
 
% get model
%==========================================================================
M         = spm_DEM_M('Lorenz');
 
% create data
%==========================================================================
 
% create innovations & add causes
%--------------------------------------------------------------------------
N         = 128;
U         = sparse(1,N);
 
% specify precisions
%--------------------------------------------------------------------------
M(1).E.nN = 32;
M(1).E.nE = 32;
M(1).E.nD = 2;
M(1).E.n  = 4;
M(1).E.d  = 2;
M(1).E.s  = 1/16;
M(1).V    = exp(1);
M(1).W    = exp(16);
DEM       = spm_DEM_generate(M,U);
 
% show data
%--------------------------------------------------------------------------
spm_figure('GetWin','Graphics'); spm_DEM_qU(DEM.pU);
 
% initialization of parameters (True values: pE =[18;-4;46.92])
%--------------------------------------------------------------------------
pE  = M(1).pE + [0 1 -1]';        % Modified prior parameters
pC  = diag(exp(-[32 2 2]));       % prior covariance
DEM.M(1).pE = pE;
DEM.M(1).pC = pC;
 
 
% Initialization of SCKS
%--------------------------------------------------------------------------
np          = length(pE);
ip          = 1:np; 
SCK         = DEM;
SCK.M(1).x  = [rand(2,1)*16;25];  % partly random initialization of the states
SCK.M(1).x  = [0.9;0.8;30];       % True state values: 
          
 
% option to specify constrains on parameters (example, not used further)
%--------------------------------------------------------------------------
cb = [0, 30;
    -10, 10;
     40, 50;];                    % [min max]
 
SCK.M(1).cb = [];                 % option to specify constrain on parameters values [min max]
SCK.M(1).ip = ip;                 % indices of model parameters to be estimated
SCK.M(2).v  = [];                 % input initial condition
SCK.M(2).V  = [];                 % input noise precision (fixed)
SCK.M(1).xP = eye(3)*1e-1^2;      % state error covariance matrix
SCK.M(1).uP = eye(1)*1e-10^2;     % input error covariance matrix
SCK.M(1).wP = pC*pC;              % parameter error covariance matrix
SCK.M(1).Q  = {speye(M(1).l,M(1).l)};     % if Q is specified then algorithm performs
                                          % estimation of measurement noise covariance 
SCK.M(1).Qf      = 'auto';        % form of estimation of measurement noise covariance 
                                  % (after online VB estimating); options:
                                  % [auto,all,min,mean]
SCK.M(1).E.nN    = 32;            % max number of iteration of SCKF-SCK algorithm
SCK.M(1).E.Itol  = 1e-3;          % convergence tolerance value for SCKF_SCK algorithm
SCK.M(1).E.RM    = [5e3 1e6];     % scaling parameter for Robbins-Monro approximation of 
                                  % parameter noise covariance [scaling parameter, max-limit]
SCK.M(1).VB.N    = 3;             % max number of VB iteration during one SCKF-SCK run
SCK.M(1).VB.Itol = 1e-2;          % convergence tolerance value for VB algorithm
SCK.M(1).VB.l    = 1 - exp(-2);   % scaling parameter for VB algorithm, 
                                  
 
% Inversion: 
%==========================================================================
LAP     = spm_LAP(DEM);
DEM     = spm_DEM(DEM);
SCK     = spm_SCK(SCK);
 
% Show estimates of states
%--------------------------------------------------------------------------
spm_figure('GetWin','spm_DEM'); spm_DEM_qU(DEM.qU,LAP.pU)
spm_figure('GetWin','spm_SCK'); spm_DEM_qU(SCK.qU,LAP.pU)
spm_figure('GetWin','spm_LAP'); spm_DEM_qU(LAP.qU,LAP.pU)
 
 
% and parameters
%--------------------------------------------------------------------------
pP1    = spm_vec(LAP.pP.P);
qP2    = spm_vec(LAP.qP.P);
qP3    = spm_vec(DEM.qP.P);
qP4    = spm_vec(SCK.qP.P);
dP1    = pE(ip)  - pP1(ip);
dP2    = qP2(ip) - pP1(ip);
dP3    = qP3(ip) - pP1(ip);
dP4    = qP4(ip) - pP1(ip);
 
 
subplot(2,2,3)
bar([pP1 qP2 qP3 qP4])
axis square
legend('True','LAP','DEM','SCK')
title('parameters (minus true)','FontSize',16)
 
subplot(2,2,4)
bar([dP1 dP2 dP3 dP4])
axis square
legend('Prior','LAP','DEM','SCK')
title('difference from true value','FontSize',16)



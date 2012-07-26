function DEM_demo_convolution_LAP
% Linear convolution revisited: A dual estimation problem
%__________________________________________________________________________
% This demonstration compares generalised filtering and a state-of-the-art 
% Bayesian smoother (SCKS) in the context of dual estimation. Note that the
% parameter estimates are smaller then the true values for generalised 
% schemes (LAP and DEM). This is largely due to the shrinkage priors and 
% optimisation of model evidence (marginal likelihood), as opposed to the
% likelihood optimised by the Square-root Cubature Kalman Smoother (SCKS).
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_convolution_LAP.m 4804 2012-07-26 13:14:18Z karl $
 
 
% create model and data
%==========================================================================
M       = spm_DEM_M('convolution model');
 
% SCK needs to be able to evaluate model equations in matrix form, i.e. for all
% cubature point at once! (in order to avoid FOR loops and speed up the process)
%--------------------------------------------------------------------------
M(1).f  = inline('[P(1,:).*x(1,:)+P(3,:).*x(2,:)+P(13,:).*v(:,:);P(2,:).*x(1,:)+P(4,:).*x(2,:)+P(14,:).*v(:,:);]','x','v','P');
M(1).g  = inline('[P(5,:).*x(1,:)+P(9,:).*x(2,:);P(6,:).*x(1,:)+P(10,:).*x(2,:);P(7,:).*x(1,:)+P(11,:).*x(2,:);P(8,:).*x(1,:)+P(12,:).*x(2,:);]','x','v','P');
M(1).V  = exp(8);                              % error precision
M(1).W  = exp(6);                              % error precision
M(2).V  = exp(8);
 
% free parameters
%--------------------------------------------------------------------------
P       = spm_vec(M(1).pE);                    % true parameters
ip      = [2 4 5 9];                           % free parameters
pE      = P;
np      = length(pE);
pE(ip)  = 0;
pC      = sparse(ip,ip,32,np,np);
M(1).pE = pE;
M(1).pC = pC;
 
% generate data
%==========================================================================
M(1).E.nE  = 32;                                % DEM-steps
M(1).E.nN  = 32;                                % DEM-steps
N          = 32;                                % length of data sequence
U          = exp(-((1:N) - 12).^2/(2.^2));      % this is the Gaussian cause
DEM        = spm_DEM_generate(M,U,{P});
DEM.M(2).V = exp(0);
 
% plot data
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); 
spm_DEM_qU(DEM.pU);
 
% Initialization of SCKS
%==========================================================================
SCK         = DEM;
SCK.M(1).cb = [];               % option to specify constrain on parameters values [min max]
SCK.M(1).ip = ip;               % indices of model parameters to be estimated
SCK.M(2).v  = 0;                % input initial condition
SCK.M(2).V  = 10;               % input noise precision (fixed)
SCK.M(1).xP = eye(2)*1e-1^2;    % state error covariance matrix
SCK.M(1).uP = eye(1)*1e-1^2;    % input error covariance matrix
SCK.M(1).wP = sparse(ip,ip,10e-10^2,np,np); % parameter error covariance matrix
SCK.M(1).pC = sparse(ip,ip,1e-3,np,np);  % covariance matrix of parameter noise
SCK.M(1).Q  = {speye(M(1).l,M(1).l)};    % if Q is specified then algorithm performs
                                         % estimation of measurement noise covariance 
SCK.M(1).Qf      = 'auto';      % form of estimation of measurement noise covariance 
                                % (after online VB estimation); options:
                                % [auto,all,min,mean]
SCK.M(1).E.nN    = 50;          % max number of iteration of SCKF-SCK algorithm
SCK.M(1).E.Itol  = 1e-3;        % convergence tolerance value for SCKF_SCK algorithm
SCK.M(1).E.RM    = [5e2 1e6];   % scaling parameter for Robbins-Monro approximation of 
                                % parameter noise covariance [scaling
                                % parameter, max-limit]
SCK.M(1).VB.N    = 10;          % max number of VB iteration during one SCKF-SCK run
SCK.M(1).VB.Itol = 1e-6;        % convergence tolerance value for VB algorithm
SCK.M(1).VB.l    = 1 - exp(-2); % scaling parameter for VB algorithm, 
                                % controls dynamics
%--------------------------------------------------------------------------
 
 
% Inversion:
%==========================================================================
DEM     = spm_DEM(DEM);
LAP     = spm_LAP(DEM);
SCK     = spm_SCK(SCK);
 
% Show estimates of states
%--------------------------------------------------------------------------
spm_figure('GetWin','spm_DEM'); spm_DEM_qU(DEM.qU,LAP.pU)
spm_figure('GetWin','spm_SCK'); spm_DEM_qU(SCK.qU,LAP.pU)
spm_figure('GetWin','spm_LAP'); spm_DEM_qU(LAP.qU,LAP.pU)
 
 
% and parameters
%--------------------------------------------------------------------------
qP1    = spm_vec(LAP.pP.P);
qP2    = spm_vec(LAP.qP.P);
qP3    = spm_vec(DEM.qP.P);
qP4    = spm_vec(SCK.qP.P);
qP1    = qP1(ip);
qP2    = qP2(ip);
qP3    = qP3(ip);
qP4    = qP4(ip);
 
subplot(2,2,4)
bar([qP1 qP2 qP3 qP4])
axis square
legend('True','LAP','DEM','SCK')
title('parameters','FontSize',16)



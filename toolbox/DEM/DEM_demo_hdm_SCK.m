function DEM_demo_hdm_SCK
% Demo for Hemodynamic deconvolution: Cross-validation of Laplace scheme
%__________________________________________________________________________
% This demonstration compares generalised filtering and SCKS in the context 
% of a nonlinear convolution model using synthetic data. Here, we look at
% estimating three of the hemodynamic parameters. This is a particularly
% difficult (almost impossible) problem, given their distance from the data
% and the conditional dependencies with the hidden states. Furthermore, 
% this is an unrealistic simulation, because we assume the data are almost 
% noiseless. The key thing to focus on is the comparative performance in
% recovering the hidden states and causes.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_hdm_SCK.m 4804 2012-07-26 13:14:18Z karl $
 
 
% create model and generate data
%==========================================================================
 
% set options
%--------------------------------------------------------------------------
clear M
M(1).E.linear = 0;                          % linear model
M(1).E.s      = 1;                          % smoothness
 
% level 1
%--------------------------------------------------------------------------
ip      = [1 2 5]';
pE      = spm_hdm_priors(1,5);              % parameters
np      = length(pE);
pC      = sparse(ip,ip,exp(-3),np,np);
pE(6)   = 0.02;
pE(7)   = 0.5;
 
M(1).n  = 4;
M(1).f  = 'spm_fx_hdm';
M(1).g  = 'spm_gx_hdm';
M(1).pE = pE;                               % prior expectation
M(1).pC = pC;                               % prior expectation
M(1).xP = exp(4);                           % prior expectation
M(1).V  = exp(8);                           % error precision
M(1).W  = exp(12);                          % error precision
M(1).ip = ip;
 
% level 2
%--------------------------------------------------------------------------
M(2).l  = 1;                                % inputs
M(2).V  = exp(0);                           % with shrinkage priors
 
% true parameters
%--------------------------------------------------------------------------
P       = M(1).pE;
P(ip)   = P(ip) - P(ip)/8;
 
 
% generate data
%==========================================================================
N         = 64;                             % length of data sequence
U         = exp(-([1:11] - 6).^2/(2.^2))/8; % this is the Gaussian cause
input     = zeros(1,N);
input([5 10 20 34 43 50]) = [1 0.8 1 0.2 .9 0.4];
U         = conv(U,input);
U         = U(1:N);
DEM       = spm_DEM_generate(M,U,{P},{8,16},{16});
 
spm_figure('GetWin','Simulated');
spm_DEM_qU(DEM.pU)
 
% Initialise SCK
%==========================================================================
DEM.M(1).E.nN = 32;
DEM.M(1).E.nN = 32;
SCK           = DEM;
 
% option to specify constrains on parameters (example, but not used here)
%-------------------------------------------------------------------------
cb(1,:) = [0.5,0.8];              % low and high bound
cb(2,:) = [0.3,0.5];              % low and high bound
cb(3,:) = [0.7,1.3];              % low and high bound
cb(4,:) = [0.29,0.35];            % low and high bound
cb(5,:) = [0.32,.4];              % low and high bound
cb(6,:) = [0.01,0.04];            % low and high bound
cb(7,:) = [0.5,0.6];              % low and high bound
 
SCK.M(1).ip = ip;  % indices of model parameters to be estimated
SCK.M(1).cb = [];  % option to specify constrain on parameters values [min max]
SCK.M(2).v  = 0;   % input initial condition
SCK.M(2).V  = 5;   % input noise precision (fixed)
SCK.M(1).xP = eye(4)*1e-1^2;      % state error covariance matrix
SCK.M(1).uP = eye(1)*1e-2^2;      % input error covariance matrix
SCK.M(1).wP = eye(np)*1e-10^2;    % parameter error covariance matrix
SCK.M(1).pC = eye(np)*1e-10^4;    % parameter error covariance matrix
SCK.M(1).f  = 'spm_fx_hdm_sck';   % state equations rewritten for matrix operations
SCK.M(1).g  = 'spm_gx_hdm_sck';   % observation equations rewritten for matrix operations
SCK.M(1).Q  = {}; 
SCK.M(1).E.nN    = 32;            % max number of iteration of SCKF-SCK algorithm
SCK.M(1).E.nD    = 4;             % time-steps
SCK.M(1).E.Itol  = 1e-3;          % convergence tolerance value for SCKF_SCK algorithm
SCK.M(1).E.RM    = [5e3 1e6];     % scaling parameter for Robbins-Monro approximation 



% % Free-energy landscape 
% %========================================================================
% spm_dem_F(DEM,ip(1));


% Inversion 
%==========================================================================
SCK     = spm_SCK(SCK);
LAP     = spm_LAP(DEM);
DEM     = spm_DEM(DEM);
 
 
% Show estimates of states
%--------------------------------------------------------------------------
spm_figure('GetWin','spm_DEM'); spm_DEM_qU(DEM.qU,LAP.pU)
spm_figure('GetWin','spm_LAP'); spm_DEM_qU(LAP.qU,LAP.pU)
spm_figure('GetWin','spm_SCK'); spm_DEM_qU(SCK.qU,LAP.pU)
 
 
% and parameters
%--------------------------------------------------------------------------
pP1    = spm_vec(LAP.M(1).pE);
qP1    = spm_vec(LAP.pP.P) - pP1;
qP2    = spm_vec(LAP.qP.P) - pP1;
qP3    = spm_vec(DEM.qP.P) - pP1;
qP4    = spm_vec(SCK.qP.P) - pP1;

 
spm_figure('GetWin','spm_DEM'); subplot(2,2,4)
bar([qP1 qP2 qP3 qP4])
axis square
legend('True','LAP','DEM','SCK')
title('parameters (minus prior)','FontSize',16)


spm_figure('GetWin','spm_SCK'); subplot(2,1,2)
t  = 1:N;
plot(t,LAP.pU.v{2},t,LAP.qU.v{2},t,DEM.qU.v{2},t,SCK.qU.v{2})
box off
legend('True','LAP','DEM','SCK')
title('Hidden causes','FontSize',16)
xlabel('time','FontSize',12)
 

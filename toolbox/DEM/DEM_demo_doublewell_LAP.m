function DEM_demo_doublewell_LAP
% The double-well revisited:
%__________________________________________________________________________
% This demonstration compares generalised filtering and a state-of-the-art 
% Bayesian smoother (SCKS) in the context of a double-well system. Here the
% Cubature filtering outperforms generalised schemes that are confounded by 
% the failure of the Laplace assumption. Note that generalised filtering
% and DEM give the same conditional estimates of states because there are 
% no free parameters or hyperparameters and the mean-field assumption 
% implcit in DEM is irrelevant.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_doublewell_LAP.m 4804 2012-07-26 13:14:18Z karl $
 
 
% Create model and data
%==========================================================================
M      = spm_DEM_M('double well');
 
% generate data (output)
%--------------------------------------------------------------------------
N         = 120;
U         = 8*sin(pi*[1:N]/16);     
M(1).E.s  = 2;
M(1).E.n  = 3;
M(1).E.nD = 4;
M(1).V    = exp(6);
M(1).W    = exp(6);
DEM       = spm_DEM_generate(M,U);
 
% show data (output)
%--------------------------------------------------------------------------
spm_figure('GetWin','Simulated'); spm_DEM_qU(DEM.pU);
 
% Initialization of SCKS
%==========================================================================
SCK         = DEM;
SCK.M(1).cb = [];          % option to specify constrain on parameters values [min max]
SCK.M(1).ip = [];          % indices of model parameters to be estimated
SCK.M(2).v  = 0;           % input initial condition
SCK.M(2).V  = 0.5;         % input noise precision (fixed)
SCK.M(1).xP = eye(1)*1e-1^2;          % state error covariance matrix
SCK.M(1).uP = eye(1)*1e0^2;           % input error covariance matrix
SCK.M(1).Q  = {speye(M(1).l,M(1).l)}; % if Q is specified then algorithm performs
                                      % estimation of measurement noise covariance
SCK.M(1).Qf      = 'auto'; % form of estimation of measurement noise covariance 
                           % (after online VB estimation); options: [auto,all,min,mean]        
SCK.M(1).E.nN    = 50;     % max number of iterations of SCKF-SCK algorithm
SCK.M(1).E.Itol  = 1e-3;   % convergence tolerance value for SCKF_SCK algorithm
SCK.M(1).VB.N    = 2;      % max number of VB iteration during one SCKF-SCK run
SCK.M(1).VB.Itol = 1e-2;   % convergence tolerance value for VB algorithm
SCK.M(1).VB.l    = 1 - exp(-1.5);     % scaling parameter for VB algorithm, 
                                      % controls dynamics
%--------------------------------------------------------------------------
 
 
% Deconvolution
%==========================================================================
SCK     = spm_SCK(SCK);
DEM     = spm_DEM(DEM);
LAP     = spm_LAP(DEM);
  
% Shoe estimates of states
%--------------------------------------------------------------------------
spm_figure('GetWin','spm_DEM'); spm_DEM_qU(DEM.qU,LAP.pU)
spm_figure('GetWin','spm_LAP'); spm_DEM_qU(LAP.qU,LAP.pU)
spm_figure('GetWin','spm_SCK'); spm_DEM_qU(SCK.qU,LAP.pU)

function DEM_demo_DEM
% Triple estimation of states, parameters and hyperparameters: This demo
% focuses on estimating both the states and parameters to furnish a
% complete system identification, given only the form of the system and its
% responses to unknown input (c.f., DEM_demo_EM, which uses known inputs)
 
% get basic convolution model
%==========================================================================
rng(1)
M       = spm_DEM_M('convolution model');
 
% free parameters
%--------------------------------------------------------------------------
P       = M(1).pE;                              % true parameters
ip      = [2,5];                                % free parameters
pE      = spm_vec(P);
np      = length(pE);
pE(ip)  = 0;
pE      = spm_unvec(pE,P);
pC      = sparse(ip,ip,exp(8),np,np);
M(1).pE = pE;
M(1).pC = pC;
 
% free hyperparameters
%--------------------------------------------------------------------------
M(1).hE = 8;
M(1).hC = 0;

M(1).gE = 0;
M(1).gC = 1;

% generate data and invert
%==========================================================================
M(1).E.nE = 32;                                % DEM-steps
N         = 32;                                % length of data sequence
U         = exp(-((1:N) - 12).^2/(2.^2));      % this is the Gaussian cause
DEM       = spm_DEM_generate(M,U,{P},{8,8},{8});
DEM       = spm_DEM(DEM);

% state estimates (inference)
%--------------------------------------------------------------------------
spm_figure('GetWin','States'); clf
spm_DEM_qU(DEM.qU,DEM.pU)

% parameter estimates (learning)
%--------------------------------------------------------------------------
spm_figure('GetWin','Parameters'); clf
spm_DEM_qP(DEM.qP,DEM.pP)

% parameter estimates (uncertainty)
%--------------------------------------------------------------------------
spm_figure('GetWin','Precisions'); clf
spm_DEM_qH(DEM.qH,DEM.pH)

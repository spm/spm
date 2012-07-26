function DEM_demo_OU
% DEM demo for linear deconvolution:  This demo considers the deconvolution
% of one of the simplest dynamical process; a random walk or Ornstein-
% Uhlenbeck process.  It shows how DEM can infer on the causes as stochastic
% innovations (c.f., Bayesian filtering) by exploiting temporal
% correlations.  Strictly speaking this is not a Ornstein-Uhlenbeck process
% because the innovations are themselves correlated and would normally be a
% Wiener process
 
% get a simple convolution model
%==========================================================================
M       = spm_DEM_M('OU');
 
% and generate data
%==========================================================================
N       = 64;                                 % length of data sequence
DEM     = spm_DEM_generate(M,N,{},{[] 8});
 
% invert model
%==========================================================================
DEM     = spm_DEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)

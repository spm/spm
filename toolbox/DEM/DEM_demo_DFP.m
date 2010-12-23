% DEM demo for linear deconvolution:  This demo considers the deconvolution
% of the responses of a single-input-multiple output input-state-output
% model (DCM) to disclose the input or causes.  It starts by demonstrating
% Variational filtering with spm_DFP; this is a stochastic filtering scheme
% that propagates particles over a changing variational energy landscape such
% that their sample density can be used to approximate the underlying
% ensemble or conditional density.  We then repeat the inversion using spm_DEM
% (i.e., under a Laplace assumption) which involves integrating the path of
% just one particle (i.e., the mode).
 
% basic deconvolution
%==========================================================================
spm_figure('GetWin','Graphics');
 
% get a simple convolution model
%==========================================================================
M      = spm_DEM_M('convolution model');
M(1).V = exp(8);
M(1).W = exp(16);

 
% and generate data
%==========================================================================
N     = 32;                                        % length of data sequence
U     = exp(-([1:N] - N/4).^2/(2*(N/32)^2));       % Gaussian cause
DEM   = spm_DEM_generate(M,U,{},{32 16});
 
% display
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.pU)
 
 
% invert model - VF
%==========================================================================
DEM  = spm_DFP(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_figure('GetWin','Graphics');
spm_DEM_qU(DEM.qU,DEM.pU)

 
% invert model - DEM
%==========================================================================
DEM  = spm_DEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)

function spm_induced_optimise_parameters
% Demo routine that optimises free parameters
%==========================================================================
%
% This exemplar routine illustrates how one can adjust or tune prior
% parameter expectations to produce desired spectral responses as specified by
% the complex eigenvalue spectrum - or a reduced form that considered a
% small number of complex values and a single (unstable) real mode.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_induced_optimise_parameters.m 4936 2012-09-18 19:47:55Z karl $
 
 
% Model specification
%==========================================================================
 
% options
%--------------------------------------------------------------------------
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = 'CMM';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
 
 
% get priors
%--------------------------------------------------------------------------
[pE pC] = spm_dcm_neural_priors({0 0 0},{},1,options.model);
[pE pC] = spm_L_priors(M.dipfit,pE,pC);
[pE pC] = spm_ssr_priors(pE,pC);
[x,f]   = spm_dcm_x_neural(pE,options.model);
 
% hidden neuronal states of interest
%--------------------------------------------------------------------------
pE.J(1:4) = [0 1 0 0];
 
 
% orders and model
%==========================================================================
nx      = length(spm_vec(x));
nu      = Ns;
u       = sparse(1,nu);
 
% create LFP model
%--------------------------------------------------------------------------
M.f     = f;
M.g     = 'spm_gx_erp';
M.x     = x;
M.n     = nx;
M.pE    = pE;
M.pC    = pC;
M.hE    = 8;
M.hC    = exp(-4);
M.m     = nu;
M.l     = Nc;
 
% solve for steady state
%--------------------------------------------------------------------------
M.x     = spm_dcm_neural_x(pE,M);
M.u     = u;
M.Hz    = 4:64;
 
 
% Target spectrum (Y) - predominantly beta and gamma
%==========================================================================
Y      = -[32 16 64 16]' + 1i*2*pi*[48 24 12 0]';
 
spm_figure('GetWin','Optimisation');
subplot(2,1,1)
[g,w]  = spm_s2csd(Y);
plot(w,g)
title('Target spectral density','FontSize',16)
xlabel('Frequency')
ylabel('Power')
 
% create generative model of the eigenspectrum (M) and invert
%--------------------------------------------------------------------------
M.IS  = 'spm_ssm2s';
for i = 1:2
    Ep   = spm_nlsi_GN(M,[],Y);
    M.pE = Ep;
end
 
% Show results with optimised parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Optimisation');
subplot(2,1,2)
[g,w]  = spm_s2csd(spm_ssm2s(Ep,M,[]));
plot(w,g)
title('predicted spectral density','FontSize',16)
xlabel('Frequency')
ylabel('Power')









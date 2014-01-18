function spm_Granger_demo
% Demo routine for induced responses
%==========================================================================
%
% This demonstration illustrates the generative or forward model used for
% time frequency responses - in other words, a biophysically plausible
% dynamic causal model for induced responses. The basic idea is to
% integrate a neural mass model to obtain expected hidden neuronal
% states produced by an unknown (parameterised) exogenous input. The states 
% are used as the expansion point for a linear perturbation analysis of
% the frequency response properties that are local in peristimulus time.
% The ensuing spectra (induced complex cross spectra) are then convolved
% with a wavelet window to generate predictions of a conventional time
% frequency (wavelet) transform. Crucially, these predictions are complex
% and can be used to characterise delays – in terms of cross covariance
% functions. Nonlinearities in the neural mass model mean that the spectral
% responses caused by random neuronal fluctuations are state dependent and
% therefore change with the expected hidden states over peristimulus time.
%
% This routine first creates a simple – two source – generative model using
% a canonical microcircuit architecture and conductance based dynamics. It
% then produces predictions of induced responses to a short and sustained
% input to the first source – as measured by two local field potential
% recordings at each source. Exactly the same model is then integrated in
% time,  using (serially correlated) random fluctuations to drive each
% source (in addition to the exogenous input). This is repeated over 16
% trials and the simulated responses are characterised in terms of a
% wavelet transform – to produce complex cross spectral  data features. 
% These are shown graphically with their analytic predictions from the 
% generative model.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_Granger_demo.m 5837 2014-01-18 18:38:07Z karl $
 
 
% Model specification
%==========================================================================
 
% number of regions
%--------------------------------------------------------------------------
N     = 1024;                                    % unmber of time bins
Nc    = 2;                                       % number of channels
Ns    = 2;                                       % number of sources
Hz    = 4:128;                                   % frequency
dt    = 4/1000;                                  % time bins
p     = 8;                                      % autoregression order
options.spatial  = 'LFP';
options.model    = 'CMM';
options.analysis = 'CSD';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
 
% extrinsic connections (forward an backward)
%--------------------------------------------------------------------------
A{1} = [0 0; 1 0];
A{2} = [0 1; 0 0];
A{3} = [0 0; 0 0];
C    = sparse(2,0);
 
 
% get priors
%--------------------------------------------------------------------------
pE    = spm_dcm_neural_priors(A,{},C,options.model);
pE    = spm_L_priors(M.dipfit,pE);
pE    = spm_ssr_priors(pE);
[x,f] = spm_dcm_x_neural(pE,options.model);

pE.J(1:4) = [0 1 0 0];

pE.a(2,:) = -1;
pE.b(1,:) = 0;
pE.c(1,:) = 0;
 
 
% orders and model
%==========================================================================
nx    = length(spm_vec(x));
 
% create forward model
%--------------------------------------------------------------------------
M.f   = f;
M.g   = 'spm_gx_erp';
M.x   = x;
M.n   = nx;
M.pE  = pE;
M.m   = Ns;
M.l   = Nc;
M.Hz  = Hz;
M.Rft = 4;


% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
M.u   = sparse(Ns,1);
 
% solve for steady state
%--------------------------------------------------------------------------
M.x   = spm_dcm_neural_x(pE,M);


% Integrate system to see response (time-series)
%==========================================================================

% Analytic spectral chararacterisation
%==========================================================================
[csd,Hz,mtf] = spm_csd_mtf(pE,M);
csd          = csd{1};
mtf          = mtf{1};
ccf          = spm_csd2ccf(csd,Hz,dt);
[mar,lag]    = spm_ccf2mar(ccf,p);

xmar.lag       = lag;
xmar.noise_cov = eye(Ns,Ns);
xmar           = spm_mar_spectra(xmar,Hz,1/dt);
for i = 1:Ns
    C(i,i) = real(sum(csd(:,i,i))/sum(xmar.P(:,i,i)));
end
xmar.noise_cov = C;
xmar           = spm_mar_spectra(xmar,Hz,1/dt);

% Numerical spectral chararacterisation
%==========================================================================

% Get spectral profile of fluctuations and noise
%--------------------------------------------------------------------------
[Gu,Gs,Gn]   = spm_csd_mtf_gu(pE,M.Hz);


% Integrate with pink noise process
%--------------------------------------------------------------------------
U.dt  = dt;
U.u   = spm_rand_power_law(Gu,Hz,dt,N);
LFP   = spm_int_L(pE,M,U);

% and measurement noise
%--------------------------------------------------------------------------
En    = spm_rand_power_law(Gn,Hz,dt,N);
Es    = spm_rand_power_law(Gs,Hz,dt,N);
E     = Es + En*ones(1,Ns);

% and estimate spectral features under a MAR model
%--------------------------------------------------------------------------
XMAR           = spm_mar(LFP + E,p);
XMAR           = spm_mar_spectra(XMAR,Hz,1/dt);
XMAR.noise_cov = C;
XMAR           = spm_mar_spectra(XMAR,Hz,1/dt);
CSD = XMAR.P;
CSP = xmar.P;


% Graphics
%==========================================================================
spm_figure('GetWin','Figure 1'); clf

% compare analytic and numerical spctral densities
%--------------------------------------------------------------------------
subplot(2,2,1)
plot((1:N)*dt,LFP)
xlabel('time')
title('LFP','FontSize',16)
axis square

subplot(2,2,2)
plot([spm_vec(xmar.lag),spm_vec(XMAR.lag)])
xlabel('expected')
ylabel('estimated')
title('auto-regression coefficients','FontSize',16)
axis square

subplot(2,2,3)
plot(Hz,abs(csd(:,1,1)),Hz,abs(CSD(:,1,1)),'--',Hz,abs(CSP(:,1,1)),':')
xlabel('frequency')
ylabel('absolute value')
title('auto-spectral density','FontSize',16)
axis square

subplot(2,2,4)
plot(Hz,abs(csd(:,2,1)),Hz,abs(CSD(:,2,1)),'--',Hz,abs(CSP(:,2,1)),':')
xlabel('frequency')
ylabel('absolute value')
title('cross-spectral density','FontSize',16)
axis square
legend('expected','estimated','expected (MAR)')

%  comparison of expected and numerical results
%==========================================================================
spm_figure('GetWin','Figure 2'); clf

% compare analytic and numerical spctral densities
%--------------------------------------------------------------------------
dtf  = xmar.dtf;
gew  = xmar.gew;
DTF  = XMAR.dtf;
GEW  = XMAR.gew;
spm_spectral_plot(Hz,mtf,'b',  'frequency','density')
spm_spectral_plot(Hz,dtf,'g',  'frequency','density')
spm_spectral_plot(Hz,DTF,'g-.','frequency','density')
spm_spectral_plot(Hz,gew,'r',  'frequency','density')
spm_spectral_plot(Hz,GEW,'r-.','frequency','density')

legend('modulation transfer function','directed transfer function','estimate','Granger causality','estimate')



 

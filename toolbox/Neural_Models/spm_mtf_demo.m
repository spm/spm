function spm_mtf_demo
% Demo routine for inverting local field potential models
%==========================================================================
%
% This demonstrates the inversion of a simple DCM for spectral activity in
% a single-source under steady-state assumptions; we use data reported
% in:
% 
% Bayesian estimation of synaptic physiology from the spectral responses of
% neural masses
% Moran, R.J.1, Stephan K.E.,  Kiebel S.J., Rombach N., O'Connor
% W.T., Murphy K.J., Reilly R.B., Friston K.J.
% 
% Abstract
% We describe a Bayesian inference scheme for quantifying the active
% physiology of neuronal ensembles using local field recordings of synaptic
% potentials. This entails the inversion of a generative neural mass model
% of steady-state spectral activity. The inversion uses Expectation
% Maximization (EM) to furnish the posterior probability of key synaptic
% parameters and the marginal likelihood of the model itself. The neural
% mass model embeds prior knowledge pertaining to both the anatomical
% [synaptic] circuitry and plausible trajectories of neuronal dynamics.
% This model comprises a population of excitatory pyramidal cells, under
% local interneuron inhibition and driving excitation from layer IV
% stellate cells. Under quasi-stationary assumptions, the model can predict
% the spectral profile of local field potentials (LFP).  This means model
% parameters can be optimised given real electrophysiological observations.
% The validity of inferences about synaptic parameters is demonstrated
% using simulated data and experimental recordings from the medial
% prefrontal cortex of control and isolation-reared Wistar rats.
% Specifically, we examined the maximum a posteriori estimates of
% parameters describing synaptic function in the two groups and tested
% predictions derived from concomitant microdialysis measures. The
% modelling of the LFP recordings revealed (i) a sensitization of
% post-synaptic excitatory responses, particularly marked in pyramidal
% cells, in the medial prefrontal cortex of socially isolated rats and (ii)
% increased neuronal adaptation.  These inferences were consistent with
% predictions derived from experimental microdialysis measures of
% extracellular glutamate levels.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mtf_demo.m 7679 2019-10-24 15:54:07Z spm $
 
 
% empirical data - sort and decimate
%--------------------------------------------------------------------------
load  'ten_minute_average_control.mat';
model = 'LFP';
y     = G_control(:);
w     = f_Control(:);
[w i] = sort(w);
y     = y(i);
 
for i = 1:64
    [d j] = min(abs(w - i));
    k(i)  = j;
end
k     = k(w(k) > 2 & w(k) < 64);
w     = w(k);                  % frequency
y     = y(k);                  % power
 
 
% specify model
%==========================================================================
 
% number of regions in coupled map lattice
%--------------------------------------------------------------------------
n     = 1;
 
% specify network (connections)
%--------------------------------------------------------------------------
A{1}  = triu(ones(n,n),1);
A{2}  = sparse(n,n);
A{3}  = sparse(n,n);
B     = {};
C     = sparse(n,1,1,n,1);
 
% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_dcm_neural_priors(A,B,C,model);

% augment with priors on spatial model
%--------------------------------------------------------------------------
[pE,pC] = spm_L_priors(n,pE,pC);

% augment with priors on endogenous inputs (neuronal) and noise
%--------------------------------------------------------------------------
[pE,pC] = spm_ssr_priors(pE,pC);

% intial states and equations of motion
%--------------------------------------------------------------------------
[x,f]  = spm_dcm_x_neural(pE,model);


% create LFP model
%--------------------------------------------------------------------------
M.dipfit.type = 'LFP';

M.IS  = 'spm_csd_mtf';
M.FS  = 'spm_fs_csd';
M.g   = 'spm_gx_erp';
M.f   = f;
M.x   = x;
M.n   = length(x);
M.pE  = pE;
M.pC  = pC;
M.hE  = 8;
M.hC  = 1/128;
M.m   = n;
M.l   = 1;
M.Hz  = w;
 
 
% inversion (in frequency space)
%==========================================================================

% data
%--------------------------------------------------------------------------
y     = spm_cond_units(y);
Y.y   = {y}; 
 
% invert
%--------------------------------------------------------------------------
Ep    = spm_nlsi_GN(M,[],Y);
 
% plot spectral density 
%--------------------------------------------------------------------------
[G w]  = spm_csd_mtf(Ep,M);
 
subplot(2,1,1)
plot(w,real(G{1}),w,y,':')
xlabel('frequency (Hz)')
xlabel('Power')
legend({'predicted','observed'})
title('Spectral inversion','FontSize',16)
axis square
grid on


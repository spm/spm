function spm_csd_demo
% Demo routine for inverting local field potential models using
% cross-spectral density summaries of steady-state dynamics
%__________________________________________________________________________
% 
% This demo illustrates the inversion of neural-mass models (Moran et al
% 2005) of steady state-responses summarised in terms of the cross-spectral
% density. These data features are extracted using a vector
% auto-regression model and transformed into frequency space for subsequent
% inversion using a biophysical neural-mass model that is parameterised in
% terms of coupling and time constants.
%
% One can generate exemplar data by integrating the neural-mass model or by
% generating data directly from the cross-spectral DCM. In this demo we 
% use the former. DCM inversion using the standard nonlinear system 
% identification scheme spm_nlsi_N (a EM-like variational scheme under the 
% Laplace assumption).
% 
% NeuroImage. 2007 Sep 1;37(3):706-20.
% A neural mass model of spectral responses in electrophysiology.Moran RJ,
% Kiebel SJ, Stephan KE, Reilly RB, Daunizeau J, Friston KJ. 
%
% Abstract:
% We present a neural mass model of steady-state membrane potentials
% measured with local field potentials or electroencephalography in the
% frequency domain. This model is an extended version of previous dynamic
% causal models for investigating event-related potentials in the
% time-domain. In this paper, we augment the previous formulation with
% parameters that mediate spike-rate adaptation and recurrent intrinsic
% inhibitory connections. We then use linear systems analysis to show how
% the model's spectral response changes with its neurophysiological
% parameters. We demonstrate that much of the interesting behaviour depends
% on the non-linearity which couples mean membrane potential to mean
% spiking rate. This non-linearity is analogous, at the population level,
% to the firing rate-input curves often used to characterize single-cell
% responses. This function depends on the model's gain and adaptation
% currents which, neurobiologically, are influenced by the activity of
% modulatory neurotransmitters. The key contribution of this paper is to
% show how neuromodulatory effects can be modelled by adding adaptation
% currents to a simple phenomenological model of EEG. Critically, we show
% that these effects are expressed in a systematic way in the spectral
% density of EEG recordings. Inversion of the model, given such
% non-invasive recordings, should allow one to quantify pharmacologically
% induced changes in adaptation currents. In short, this work establishes a
% forward or generative model of electrophysiological recordings for
% psychopharmacological studies.
% 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd_demo.m 6236 2014-10-12 10:03:44Z karl $
 
clear global
 
% specify model
%==========================================================================
 
% number of sources and LFP channels (usually the same)
%--------------------------------------------------------------------------
n     = 2; % number of sources
nc    = 2; % number of channels
 
% specify network (connections)
%--------------------------------------------------------------------------
A{1}  = tril(ones(n,n),-1);                % a forward connection
A{2}  = triu(ones(n,n),+1);                % a backward connection
A{3}  = sparse(n,n);                       % lateral connections
B     = {};                                % trial-specific modulation
C     = speye(n,n);                        % sources receiving innovations
 
% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_lfp_priors(A,B,C);           % neuronal priors
[pE,pC] = spm_ssr_priors(pE,pC);           % spectral priors
[pE,pC] = spm_L_priors(n,pE,pC);           % spatial  priors

% Suppress channel noise
%--------------------------------------------------------------------------
pE.b  = pE.b - 16;
pE.c  = pE.c - 16;
 
% create LFP model
%--------------------------------------------------------------------------
M.dipfit.type = 'LFP';

M.IS = 'spm_csd_mtf';
M.FS = 'spm_fs_csd';
M.g  = 'spm_gx_erp';
M.f  = 'spm_fx_lfp';
M.x  = sparse(n,13);
M.n  = n*13;
M.pE = pE;
M.pC = pC;
M.m  = n;
M.l  = nc;
M.Hz = (1:64)';
 
 
% simulate spectral data directly
%==========================================================================
P           = pE;
P.A{1}(2,1) = 1/2;                          % strong forward connections
CSD         = spm_csd_mtf(P,M);
CSD         = CSD{1};
 
% or generate data and use the sample CSD
%==========================================================================
 
% Integrate with pink noise process
%--------------------------------------------------------------------------
N    = 512;
U.dt = 8/1000;
U.u  = randn(N,M.m)/16;
U.u  = sqrt(spm_Q(1/16,N))*U.u;
LFP  = spm_int_L(P,M,U);
 
% and estimate spectral features under a MAR model
%--------------------------------------------------------------------------
try
    mar = spm_mar(LFP,8);
catch
    warndlg('please include spectral toolbax in Matlab path')
end
mar  = spm_mar_spectra(mar,M.Hz,1/U.dt);


spm_figure('GetWin','Figure 1'); clf

subplot(2,1,1)
plot((1:N)*U.dt,LFP)
xlabel('time')
title('LFP')
 
subplot(2,1,2)
plot(M.Hz,real(CSD(:,1,1)),M.Hz,real(CSD(:,1,2)),':')
xlabel('frequency')
title('[cross]-spectral density')
axis square

 
% inversion (in frequency space)
%==========================================================================
 
% data and confounds
%--------------------------------------------------------------------------
Y.y   = {CSD};
 
% invert
%--------------------------------------------------------------------------
Ep    = spm_nlsi_GN(M,[],Y);
 

spm_figure('GetWin','Figure 2'); clf
 
% plot spectral density
%==========================================================================
[G w] = spm_csd_mtf(Ep,M);
 
% plot
%--------------------------------------------------------------------------
g = G{1};
y = Y.y{1};
for i = 1:nc
    for j = 1:nc
        
        subplot(3,2,(i - 1)*nc + j)
        plot(w,real(g(:,i,j)),w,real(y(:,i,j)),':')
        title(sprintf('cross-spectral density %d,%d',i,j))
        xlabel('Power')
        axis square
        
        try axis(a), catch, a = axis; end
 
    end
end
legend({'predicted','observed'})
 
% plot parameters and estimates
%--------------------------------------------------------------------------
subplot(3,2,5)
bar(exp(spm_vec(P)))
title('true parameters')
 
subplot(3,2,6)
bar(exp(spm_vec(Ep)))
title('conditional expectation')

function spm_nested_oscillations_demo
% Demo routine for neural mass models of nested oscillations
%==========================================================================
% 
% This demo simply illustrates nested oscillations in a three-subpopulation
% source that are caused by nonlinear interactions between voltage and
% conductance. Put simply, a slow sinusoidal drive elicits periods of bursting
% to produce phase-amplitude coupling in the ensuing dynamics.  We look at
% this using both neural-mass and mean-field models.  See Marreiros et al:
% 
% Population dynamics under the Laplace assumption.
% 
% A Marreiros, J Daunizeau, S Kiebel, L Harrison & Karl Friston
% 
% Abstract
% In this paper, we describe a generic approach to modelling dynamics in
% neuronal populations.  This approach retains a full density on the states
% of neuronal populations but resolves the problem of solving
% high-dimensional problems by re-formulating density dynamics in terms of
% ordinary differential equations on the sufficient statistics of the
% densities considered.  The particular form for the population density we
% adopt is a Gaussian density (c.f., a Laplace assumption). This means
% population dynamics are described completely by equations governing the
% evolution of the population’s mean and covariance.  We derive these
% equations from the Fokker-Planck formalism and illustrate their
% application to a reasonably simple conductance-based model of neuronal
% exchanges.  One interesting aspect of this formulation is that we can
% uncouple the mean and covariance to furnish a neural-mass model, which
% rests only on the populations mean.  This enables to compare equivalent
% mean-field and neural-mass models of the same populations and evaluate,
% quantitatively, the contribution of population variance to the expected
% dynamics.  The mean-field model presented here will form the basis of a
% dynamic causal model of observed electromagnetic signals in future work.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nested_oscillations_demo.m 5934 2014-03-28 15:03:00Z karl $
 
 
% number of regions in coupled map lattice
%--------------------------------------------------------------------------
n     = 1;
 
% extrinsic network connections
%--------------------------------------------------------------------------
if n > 1
A{1}  = diag(ones(n - 1,1),-1);
else
    A{1} = 0;
end
A{2}  = A{1}';
A{3}  = sparse(n,n);
B     = {};
C     = sparse(1,1,1,n,1);
 
 
% get connectivity and other priors
%--------------------------------------------------------------------------
[pE,pC] = spm_nmm_priors(A,B,C);           % neuronal priors
[pE,pC] = spm_L_priors(n,pE,pC);           % spatial  priors
 
% initialise states and models (mean-field and neural mass)
%--------------------------------------------------------------------------
[x MF] = spm_x_mfm(pE);
[x NM] = spm_x_nmm(pE);
 
 
% create exogenous inputs for responses to transient and sustained input
%==========================================================================
W     = 4;                          % frequency of exogenous input (theta)
dt    = 4;
t     = (1:dt:1024)';
U.dt  = dt/1000;
U.u   = 4*exp(sin(2*pi*W*t/1000));
 
% responses to different inputs - spikes
%==========================================================================
p     = 3;                                               % pyramidal cells
 
% Integrate systems and plot
%----------------------------------------------------------------------
MFM   = spm_int_L(pE,MF,U);
NMM   = spm_int_L(pE,NM,U);

spm_figure('GetWin','Nested Oscillations');
 
subplot(2,1,1)
plot(t,MFM(:,1:3),t,U.u,':')
xlabel('time (ms)')
ylabel('depolarisation (mV)')
title('Mean-field model')
 
subplot(2,1,2)
plot(t,NMM(:,1:3),t,U.u,':')
xlabel('time (ms)')
ylabel('depolarisation (mV)')
title('Neural-mass model')
 
legend({'spiny cells', 'inhibitory cells', 'pyramidal cells', 'input'})

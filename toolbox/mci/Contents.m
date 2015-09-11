%
% MONTE CARLO INFERENCE (MCI) toolbox
%
% INFERENCE FOR SINGLE DATA SET
%
% spm_mci_lgv.m     Langevin Monte Carlo (LMC) otherwise known as Simplified 
%                   Manifold Metropolis Adjusted Langevin Algorithm 
%                   (Simplified MMALA). 
%
% spm_mci_ais.m     Annealed Importance Sampling (AIS)
%
% spm_mci_pop.m     Adaptive Monte Carlo (AMC): Metropolis-Hastings with 
%                   proposals tuned using Robbins-Monro. Also allows for 
%                   multiple chains and thermodynamic integration
%
% spm_mci_post.m    Generic wrapper for single subject Bayesian inference.
%                   Implemented using LMC, AIS, AMC or Variational Laplace (VL) 
%
% mci_demo_approach.m   Approach to limit example
% mci_demo_discount.m   Temporal discounting model
% mci_demo_growth.m     Preece-Baines growth model
% mci_demo_lds.m        Linear dynamical system
% mci_demo_linear.m     Linear regression
% mci_demo_linsqr.m     (Squared) Linear regression with local minima
% mci_demo_logistic.m   Logistic regression
% mci_demo_nmm.m        Neural mass models
% mci_demo_phase.m      Fully connected phase coupling models
% mci_demo_rphase.m     Phase coupling models with specific connectivity
% mci_demo_ramsay.m     Nonlinear oscillator with local minima
% 
% When setting up inference for a new dynamical model, use spm_mci_check(M) to 
% see that model M has required fields.
%
% INFERENCE FOR GROUP DATA
%
% spm_mci_mfx.m             Mixed effects inference for nonlinear systems
% spm_mci_mfx_dynamic.m     Mixed effects inference for dynamical systems
%  
% mci_demo_rfx_linear.m     Random effects linear regression and comparison
%                           with parametric Empirical Bayes
% mci_demo_rfx_logistic.m   Random effects logistic regression
% mci_demo_rfx_nmm.m        Random effects neural mass models
% mci_demo_rfx_rphase.m     Random effects phase coupling
% mci_demo_mfx_lds.m        Mixed effects linear dynamical systems
%
% INTEGRATION
%
% spm_mci_fwd.m             Integrate dynamics and apply observation model.
% spm_mci_sens.m            Forward Sensitivity analysis
% spm_mci_adjoint.m         Adjoint Sensitivity analysis
%
% mci_compare_forward.m     Compare integration speed of various methods
% mci_compare_gradients.m   Compare accuracy of gradient estimation
% mci_compare_sensitivities Compare speed of sensitivity estimation
%
% The sensitivity matrices are more efficiently computed if you have 
% installed the four major components of the Sundials package (CVODE,
% CVODES,IDA,IDAS) from http://computation.llnl.gov/casc/sundials/
%
% DIAGNOSTICS
%
% spm_mci_diag.m            Trace plots, energy trajectory
% spm_mci_ess.m             Effective sample size for a Markov chain
% spm_mci_stat.m            Test for stationarity
% spm_mci_quantiles         Histograms and quantiles from samples
% 
% REFERENCES:
%
% W.Penny, M Klein-Flugge and B Sengupta (2015) Mixed-Effects Langevin 
% Monte Carlo, Submitted, 2015.
%
% W.Penny and B Sengupta (2015) Annealed Importance Sampling for Neural 
% Mass Models, Submitted, 2015.
%
% B. Sengupta, K. Friston and W. Penny (2015) Gradient-based MCMC samplers
% for dynamic causal modelling. Neuroimage. 
%
% B. Sengupta, K. Friston and W. Penny (2015) Gradient-free MCMC samplers
% for dynamic causal modelling. Neuroimage, 112, 375-381.
%
% B. Sengupta, K. Friston and W. Penny (2014) Efficient Gradient
% Computation for Dynamical Models. Neuroimage,98, 521-527. 
%_________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: Contents.m 6548 2015-09-11 12:39:47Z will $
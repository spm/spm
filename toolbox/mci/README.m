% MONTE CARLO INFERENCE (MCI) toolbox
%
% CORE FUNCTIONS:
%
% spm_mci_mfx.m     Mixed effects models for dynamical systems
%  
% spm_mci_lgv.m     Langevin Monte Carlo - Simplified Manifold Metropolis 
%                   Adjusted Langevin Algorithm (Simplified MMALA). 
%
% spm_mci_pop.m     Metropolis-Hastings (MH) with proposals tuned 
%                   using Robbins-Monro. Also allows for 
%                   multiple chains and thermodynamic integration
%
% spm_mci_post.m    Generic wrapper for single subject Bayesian inference.
%                   Implemented using MH, Variational Laplace or Langevin
%
% spm_mci_fwd.m     Integrate dynamics and apply observation model.
%
% spm_mci_sens.m    Forward Sensitivity analysis
%
% spm_mci_adjoint.m Adjoint Sensitivity analysis
%
% The sensitivity matrices are more efficiently computed if you have 
% installed the four major components of the Sundials package (CVODE,
% CVODES,IDA,IDAS) from http://computation.llnl.gov/casc/sundials/
%
% The data structures also allow for inference using Variational Laplace
% (VL). See eg spm_mci_post.m. When setting up inference for a new
% model, use spm_mci_chk(M) to see that the model M has required fields.
%
% DEMOS:
%
% fitzhugh          Fitzhugh-Nagumo model
% lds               Linear Dynamical System
% linear            Linear Model
% nmm               Neural Mass Model
% phase             Weakly Coupled Oscillator
%
% REFERENCES:
%
% B Sengupta, G Ziegler, G Ridgway, K Friston, J Ashburner and W.Penny.
% Mixed Effects Models of Dynamical Systems, Submitted, 2014.
%
% B. Sengupta, K. Friston and W. Penny (2014) Efficient Gradient
% Computation for Dynamical Models. Neuroimage,98, 521-527. 
%_________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: README.m 6275 2014-12-01 08:41:18Z will $
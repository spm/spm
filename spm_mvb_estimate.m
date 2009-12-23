function [MVB] = spm_mvb_estimate(MVB)
% [re]esimates a multivariate Bayes model (Bayesian decoding of a contrast)
% FORMAT [MVB] = spm_mvb_estimate(MVB)
%
% Sets up, evaluates and saves an MVB structure:
%
% MVB.X                   % subspace of design matrix
% MVB.Y                   % multivariate response
% MVB.X0                  % null space of design
% MVB.V                   % serial correlation in response
% MVB.XYZ                 % location of voxels (mm)
% MVB.VOX                 % voxel scaling
% MVB.Ni                  % number of greedy search steps
% MVB.sg                  % size of reedy search split
%
% where MVB.M will contain the following fields:
%
%                F: log-evidence [F(0), F(1),...]
%                G: covariance partition indices
%                h: covariance hyperparameters
%                U: ordered patterns
%               qE: conditional expectation of voxel weights
%               qC: conditional variance of voxel weights
%               Cp: prior covariance (ordered  pattern space)
%               cp: prior covariance (original pattern space)
%
%--------------------------------------------------------------------------
% This routine uses a multivariate Bayesian (MVB) scheme to decode or
% recognise brain states from neuroimages. It resolves the ill-posed
% many-to-one mapping, from voxel values or data features to a target
% variable, using a parametric empirical or hierarchical Bayesian model.
% This model is inverted using standard variational techniques, in this
% case expectation maximisation, to furnish the model evidence and the
% conditional density of the model's parameters. This allows one to compare
% different models or hypotheses about the mapping from functional or
% structural anatomy to perceptual and behavioural consequences (or their
% deficits). The aim of MVB is not to predict (because the outcomes are
% known) but to enable inference on different models of structure-function
% mappings; such as distributed and sparse representations. This allows one
% to optimise the model itself and produce predictions that outperform
% standard pattern classification approaches, like support vector machines.
% Technically, the model inversion and inference uses the same empirical
% Bayesian procedures developed for ill-posed inverse problems (e.g.,
% source reconstruction in EEG).
%
% CAUTION: MVB should not be used to establish a significant mapping
% between brain states and some classification or contrast vector. Its use
% is limited to comparison of different models under the assumption
% (hyperprior) that this mapping exists. To ensure the mapping exists, use
% CVA or compute the randomisation p-value (see spm_mvb_p)
%
% See: spm_mvb and
%
% Bayesian decoding of brain images.
% Friston K, Chu C, Mourao-Miranda J, Hulme O, Rees G, Penny W, Ashburner J.
% Neuroimage. 2008 Jan 1;39(1):181-205
% 
% Multiple sparse priors for the M/EEG inverse problem.
% Friston K, Harrison L, Daunizeau J, Kiebel S, Phillips C, Trujillo-Barreto 
% N, Henson R, Flandin G, Mattout J.
% Neuroimage. 2008 Feb 1;39(3):1104-20.
% 
% Characterizing dynamic brain responses with fMRI: a multivariate approach.
% Friston KJ, Frith CD, Frackowiak RS, Turner R.
% Neuroimage. 1995 Jun;2(2):166-72.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_estimate.m 3656 2009-12-23 20:17:30Z karl $

%-Figure
%--------------------------------------------------------------------------
spm_figure('GetWin','MVB');

%-Get model[s]
%--------------------------------------------------------------------------
str    = {'compact','sparse','smooth','support'};
Ip     = spm_input('model (spatial prior)','!+1','m',str);
priors = str{Ip};

% invert
%==========================================================================
VOX        = diag(MVB.VOX);
U          = spm_mvb_U(MVB.Y,priors,MVB.X0,MVB.XYZ,VOX(1:end - 1));
M          = spm_mvb(MVB.X,MVB.Y,MVB.X0,U,MVB.V,MVB.Ni,MVB.sg);
M.priors   = priors;
 
% assemble results
%--------------------------------------------------------------------------
MVB.M      = M;                           % MVB model (see below)
MVB.priors = priors;                      % model (spatial prior)
 
% display
%==========================================================================
spm_mvb_display(MVB)

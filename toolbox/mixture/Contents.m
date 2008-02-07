% SPM Mixture Modelling Toolbox 
%
% Bayesian Multivariate mixture modelling [1,2]: spm_mix.m
% 
% Robust General Linear Model [3]: spm_rglm.m
%
% Kmeans clustering: spm_kmeans.m, spm_kmeans1.m
%
% A number of routines are based on NETLAB functions 
% (see http://www.ncrg.aston.ac.uk/netlab/), though NETLAB is not
% required in your search path.
% 
% References:
%
% [1] H. Attias (2000) A Variational Bayesian framework for Graphical Models, 
%     NIPS 12, 209-215, MIT press, 2000.
%
% [2] W.D. Penny (2001) Variational Bayes for d-dimensional Gaussian mixture models 
%     Wellcome Department of Imaging Neuroscience, University College London.
%
% [3] W.D. Penny and J. Kilner (2007) Robust Bayesian General Linear
%     Models. Neuroimage.
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: Contents.m 1143 2008-02-07 19:33:33Z spm $
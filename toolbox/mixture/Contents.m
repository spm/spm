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
%__________________________________________________________________________
% Copyright (C) 2007-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: Contents.m 5962 2014-04-17 12:47:43Z spm $
function [pE,pC,x] = spm_dcm_nvc_priors(DCM)
% Priors for a multimodal DCM for fMRI and M/EEG
% FORMAT [pE,pC,x] = spm_dcm_nvc_priors(DCM)
%
% Input:
% -------------------------------------------------------------------------
% DCM      - multimodal DCM (see spm_dcm_nvc_specify.m)
%
% Evaluates:
% -------------------------------------------------------------------------
% pE.H     - prior expectations (hemodynamic)
% pC.H     - prior covariances  (hemodynamic)
% pE.J     - prior expectations (neurovascular coupling)
% pC.J     - prior covariances  (neurovascular coupling)
% x        - prior (initial) states
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Humman Neuroimaging

% Amirhossein Jafarian 
% $Id: spm_dcm_nvc_priors.m 7713 2019-11-25 16:00:34Z spm $

% Number of regions & model space (see spm_dcm_nvc.m for model details)
%--------------------------------------------------------------------------
n       = DCM.n;
model   = DCM.model;

% Haemodynamic priors for HDM model
%--------------------------------------------------------------------------
pE.H.transit = sparse(n,1);  pC.H.transit = sparse(n,1) + 1/256;
pE.H.decay   = sparse(n,1);  pC.H.decay   = sparse(n,1) + 1/256;
pE.H.epsilon = sparse(1,1);  pC.H.epsilon = sparse(1,1) + 1/256;
x            = sparse(n,4);

% Priors for neurovascular parameters
%--------------------------------------------------------------------------
if (strcmp(model(1), 'pre')&& strcmp(model(2), 'd') && strcmp(model(3),'int'))
    pE.J =  sparse(4,n);
    pC.J =  sparse(4,n)   + 1/16 *repmat(DCM.N',1,n);
elseif (strcmp(model(1), 'pre') && strcmp(model(2), 's') && strcmp(model(3),'int'))
    pE.J =  sparse(4,1);
    pC.J =  sparse(4,1) + 1/16*(DCM.N'); 
elseif (strcmp(model(1), 'pre')&& strcmp(model(2), 'd') && strcmp(model(3),'ext'))
    pE.J =  sparse(4,n);
    pC.J =  sparse(4,n) + 1/16*repmat(DCM.N',1,n);
elseif (strcmp(model(1), 'pre') && strcmp(model(2), 's')  &&  strcmp(model(3),'ext'))
    pE.J =  sparse(4,1);
    pC.J =  sparse(4,1) + 1/16*(DCM.N');    
elseif (strcmp(model(1), 'de')&& strcmp(model(2), 's') && strcmp(model(3),'int'))
    pE.J =  sparse(4,2);
    pC.J =  sparse(4,2) + 1/16*repmat(DCM.N',1,2);
elseif (strcmp(model(1), 'de') && strcmp(model(2), 's')&& strcmp(model(3),'ext'))
    pE.J =  sparse(4,3);
    pC.J =  sparse(4,3) + 1/16*repmat(DCM.N',1,3);   
elseif (strcmp(model(1), 'de')&& strcmp(model(2), 'd') && strcmp(model(3),'int'))
    pE.J =  zeros(4,2,n);
    pC.J =  zeros(4,2,n) + 1/16*reshape(repmat(DCM.N',2*n,1), [4 2 n]);
elseif (strcmp(model(1), 'de') && strcmp(model(2), 'd')&& strcmp(model(3),'ext'))
    pE.J =  zeros(4,3,n);
    pC.J =  zeros(4,3,n) + 1/16*reshape(repmat(DCM.N',3*n,1), [4 3 n]);   
elseif (strcmp(model(1), 'post') && strcmp(model(2), 's') &&  strcmp(model(3),'na'))
    pE.J =  sparse(4,1);
    pC.J =  sparse(4,1) + 1/16*(DCM.N'); 
elseif (strcmp(model(1), 'post') && strcmp(model(2), 'd') && strcmp(model(3),'na'))
    pE.J =  sparse(4,n);
    pC.J =  sparse(4,n) + 1/16*repmat(DCM.N',1,n);
end

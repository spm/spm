function [pE,pC,x] = spm_dcm_mm_priors(DCM)
% Priors for a multimodal DCM for fMRI and M/EEG
% FORMAT:[pE,pC,x] = spm_dcm_mm_priors(DCM)
%
%----------------------------------Input----------------------------------- 
% DCM 
%
%--------------------------------- Output----------------------------------
%    pE.H     - prior expectations (hemodynamic)
%    pC.H     - prior covariances  (hemodynamic)
%    pE.J     - prior expectations (neurovascular coupling)
%    pC.J     - prior covariances  (neurovascular coupling)
%    x        - prior (initial) states
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Humman Neuroimaging

% Amirhossein Jafarian 
% $Id: spm_dcm_mm_priors.m 7705 2019-11-22 15:06:38Z spm $


% Number of regions & Model space 
%--------------------------------------------------------------------------
n       = DCM.n;
Model   = DCM.model ;

% Haemodynamic prior for HDM model
%--------------------------------------------------------------------------
bE.transit = sparse(n,1);  bC.transit = sparse(n,1) + 1/256;
bE.decay   = sparse(n,1);  bC.decay   = sparse(n,1) + 1/256;
bE.epsilon = sparse(1,1);  bC.epsilon = sparse(1,1) + 1/256;
x          = sparse(n,4);
pE.H       = bE;
pC.H       = bC;

% Prior for neurovascular parameters
%--------------------------------------------------------------------------
if (strcmp(Model(1), 'pre')&& strcmp(Model(2), 'd') && strcmp(Model(3),'int'))
    pE.J =  sparse(4,n);
    pC.J =  sparse(4,n)   + 1/16 *repmat(DCM.N',1,n);
elseif (strcmp(Model(1), 'pre') && strcmp(Model(2), 's') && strcmp(Model(3),'int'))
    pE.J =  sparse(4,1);
    pC.J =  sparse(4,1) + 1/16*(DCM.N'); 
elseif (strcmp(Model(1), 'pre')&& strcmp(Model(2), 'd') && strcmp(Model(3),'ext'))
    pE.J =  sparse(4,n);
    pC.J =  sparse(4,n) + 1/16*repmat(DCM.N',1,n);
elseif (strcmp(Model(1), 'pre') && strcmp(Model(2), 's')  &&  strcmp(Model(3),'ext'))
    pE.J =  sparse(4,1);
    pC.J =  sparse(4,1) + 1/16*(DCM.N');    
elseif (strcmp(Model(1), 'de')&& strcmp(Model(2), 's') && strcmp(Model(3),'int'))
    pE.J =  sparse(4,2);
    pC.J =  sparse(4,2) + 1/16*repmat(DCM.N',1,2);
elseif (strcmp(Model(1), 'de') && strcmp(Model(2), 's')&& strcmp(Model(3),'ext'))
    pE.J =  sparse(4,3);
    pC.J =  sparse(4,3) + 1/16*repmat(DCM.N',1,3);   
elseif (strcmp(Model(1), 'de')&& strcmp(Model(2), 'd') && strcmp(Model(3),'int'))
    pE.J =  zeros(4,2,n);
    pC.J =  zeros(4,2,n) + 1/16*reshape(repmat(DCM.N',2*n,1), [4 2 n]);
elseif (strcmp(Model(1), 'de') && strcmp(Model(2), 'd')&& strcmp(Model(3),'ext'))
    pE.J =  zeros(4,3,n);
    pC.J =  zeros(4,3,n) + 1/16*reshape(repmat(DCM.N',3*n,1), [4 3 n]);   
elseif (strcmp(Model(1), 'post') && strcmp(Model(2), 's') &&  strcmp(Model(3),'na'))
    pE.J =  sparse(4,1);
    pC.J =  sparse(4,1) + 1/16*(DCM.N'); 
elseif (strcmp(Model(1), 'post') && strcmp(Model(2), 'd') && strcmp(Model(3),'na'))
    pE.J =  sparse(4,n);
    pC.J =  sparse(4,n) + 1/16*repmat(DCM.N',1,n);
end

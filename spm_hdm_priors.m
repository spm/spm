function [pE,pC] = spm_hdm_priors(m)
% returns priors for a hemodynamic dynmaic causal model
% FORMAT [pE,pC] = spm_hdm_priors(m)
% m.. - number of inputs
%
% pE  - prior expectations
% pC  - prior covariances
%___________________________________________________________________________
% %W% Karl Friston %E%

% append input efficacy priors to (5) biophysical parameters
%---------------------------------------------------------------------------
pE    = [0.65      0.41      0.98      0.32      0.34  ];
pC    = [0.0150    0.0020    0.0568    0.0013    0.0024];

pE    = [pE(:); zeros(m,1)];
pC    = blkdiag(diag(pC),eye(m)*16);
return

% sample covariances from Friston et al (2000)
%---------------------------------------------------------------------------
pC    = [   0.0150    0.0052    0.0283    0.0002   -0.0027
	    0.0052    0.0020    0.0104    0.0004   -0.0013
	    0.0283    0.0104    0.0568    0.0010   -0.0069
	    0.0002    0.0004    0.0010    0.0013   -0.0010
	   -0.0027   -0.0013   -0.0069   -0.0010    0.0024];

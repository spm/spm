function [pE,pC,qE,qC] = spm_dcm_priors(A,B,C,varargin)
% Returns the priors for a hemodynamic dynamic causal model.
% FORMAT:
%    for bi-linear DCM: [pE,pC,qE,qC] = spm_dcm_priors(A,B,C) 
%    for nonlinear DCM: [pE,pC,qE,qC] = spm_dcm_priors(A,B,C,D) 
% INPUT:
%    A,B,C,D - constraints on connections (1 - present, 0 - absent)
%
% OUTPUT:
%    pE     - prior expectations (connections and hemodynamic)
%    pC     - prior covariances  (connections and hemodynamic)
%    qE     - prior expectations (hemodynamic)
%    qC     - prior covariances  (hemodynamic)
% 
% NB: order of parameter vector is: 
% bilinear neural params -> hemodynamic params -> nonlinear neural params
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_priors.m 3705 2010-02-01 20:51:28Z karl $
 


% number of regions
%--------------------------------------------------------------------------
n     = length(A);       % number of regions
q     = 1/32;            % prior covariance, if A(i,j) = 1
b     = 1;               % global decay rate (Hz)

% varargin (D for nonlinear coupling)
%--------------------------------------------------------------------------
if nargin > 3, D = varargin{1}; else, D = zeros(n,n,0); end
 
% prior covariances
%--------------------------------------------------------------------------
A     = A - diag(diag(A));
pC    = diag([1/16; A(:)*q; B(:)*q; C(:) D(:)*q]);

% prior expectations 
%--------------------------------------------------------------------------
pE.b  =  log(b);
pE.A  = -speye(n,n);
pE.B  =  B*0;
pE.C  =  C*0;
pE.D  =  D*0;


% and combine with hemodynamic priors
%==========================================================================
% P(1) - signal decay     - d(ds/dt)/ds)  half-life = log(2)/P(1) ~ 1sec
% P(2) - auto-regulation  - d(ds/dt)/df)  2*pi*sqrt(1/P(2)) ~ 10 sec
% P(3) - transit time               (t0)  ~ 1 sec
% P(4) - exponent for Fout(v)    (alpha)  c.f. Grubb's exponent (~ 0.38)
% P(5) - resting oxygen extraction  (E0)  ~ range 20 - 50%
% P(6) - ratio: intra- to extravascular components of gradient echo signal:
%        epsilon (prior mean = 1, log-normally distributed scaling factor)

[qE,qC] = spm_hdm_priors(0);

% combine connectivity and hemodynamic priors
%==========================================================================
qC      = kron(qC,   eye(n,n));
qE      = kron(qE', ones(n,1));
pE.H    = qE;
pC      = blkdiag(pC,qC);

return

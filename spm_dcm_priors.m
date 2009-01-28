function [pE,pC,qE,qC] = spm_dcm_priors(A,B,C,varargin)
% Returns the priors for a hemodynamic dynamic causal model.
% FORMAT:
%    for bilinear DCM:  [pE,pC,qE,qC] = spm_dcm_priors(A,B,C) 
%    for nonlinear DCM: [pE,pC,qE,qC] = spm_dcm_priors(A,B,C,D) 
% INPUT:
%    A,B,C,D - constraints on connections (1 - present, 0 - absent)
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
% $Id: spm_dcm_priors.m 2661 2009-01-28 20:21:42Z karl $
 

% nonlinear DCM?
%--------------------------------------------------------------------------
if nargin > 3
    nlDCM = 1;
    D     = varargin{1};
else
    nlDCM = 0;
end

% number of regions
%--------------------------------------------------------------------------
n     = length(A);
 

% CONNECTIVITY PRIORS
%==========================================================================
% covariances            pC = 1/((n - 1)*q)  - if A == 1
%                        pC = 0                if A == 0
%
% where, for aij = aji = 1/(n - 1)  => max(eig(J(0))) = 0
% and q is the Chi-squared threshold with n*(n - 1) df
 
% log(2)/b = half-life {b = self inhibition} - log-normal distribution
%--------------------------------------------------------------------------
p     = 1e-3;
b     = 1;
s     = spm_invNcdf(1 - p);
q     = spm_invXcdf(1 - p,n*(n - 1));
q     = n/((n - 1)*q);
 
% intrinsic connections A {additional priors from eigenvalues}
%--------------------------------------------------------------------------
A     = A - diag(diag(A));
pC    = diag([1/16; A(:)*q; B(:)*q; C(:)]);
if ~nlDCM
    % bilinear DCM
    pC_D = [];
else
    % nonlinear DCM
    pC_D  = diag(D(:)*1);
end

% define prior expectations and combine with hemodynamic priors
%--------------------------------------------------------------------------
A     = -speye(n,n);
B     = B*0;
C     = C*0;
pE    = [log(b); A(:); B(:); C(:)];
if ~nlDCM
    % bilinear DCM
    pE_D = [];
else
    % nonlinear DCM
    D     = D*0;
    pE_D  = [D(:)];
end


% HEMODYNAMIC PRIORS
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
qC      = kron(qC,eye(n,n));
qE      = kron(qE,ones(n,1));
pE      = [pE; qE; pE_D];
pC      = blkdiag(pC,qC,pC_D);


return

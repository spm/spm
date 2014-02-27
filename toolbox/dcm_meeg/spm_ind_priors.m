function [pE,gE,pC,gC] = spm_ind_priors(A,B,C,Nm,Nf)
% prior moments for a neural-mass model of induced responses
% FORMAT [pE,gE,pC,gC] = spm_ind_priors(A,B,C,dipfit,Nu,Nf)
% A{2},B{m},C  - binary constraints on extrinsic connections
% Nm           - number of frequency modes used
% Nf           - number of frequency modes explained
%
% pE - prior expectation - f(x,u,P,M)
% gE - prior expectation - g(x,u,G,M)
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A    - trial-invariant
%    pE.B{m} - trial-dependent
%    pE.C    - stimulus-stimulus dependent
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R    - onset and dispersion
%
% pC - prior covariances: cov(spm_vec(pE))
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ind_priors.m 5900 2014-02-27 21:54:51Z karl $
 
% orders
%--------------------------------------------------------------------------
if nargin < 5; Nf = Nm; end
n    = size(C,1);                                 % number of sources
nu   = size(C,2);                                 % number of inputs
D    = sparse(1:Nm,1:Nm,1,Nf,Nf);


% set extrinsic connectivity - linear and nonlinear (cross-frequency)
%--------------------------------------------------------------------------
E.A  = kron(speye(Nf,Nf), -speye(n,n));
V.A  = kron(D,A{1}) + kron(D*(1 - speye(Nf,Nf))*D,A{2});
 
% input-dependent
%--------------------------------------------------------------------------
for i = 1:length(B)
    E.B{i} = sparse(n*Nf,n*Nf);
    V.B{i} = kron(D*ones(Nf,Nf)*D,B{i}) & V.A;
end
 
% exogenous inputs
%--------------------------------------------------------------------------
E.C  = kron(D*ones(Nf,1),C - C);
V.C  = kron(D*ones(Nf,1),C);
 
% set stimulus parameters: magnitude, onset and dispersion
%--------------------------------------------------------------------------
E.R  = kron(ones(nu,1),[0 0]);
V.R  = kron(ones(nu,1),[1/16 1/16]);
 
% prior moments
%--------------------------------------------------------------------------
pE   = E;
gE   = [];
pC   = diag(sparse(spm_vec(V)));
gC   = [];

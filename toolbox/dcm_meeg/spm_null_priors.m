function [pE,pC] = spm_null_priors(A,B,C)
% prior moments for null (Jacobian) model
% FORMAT [pE,pC] = spm_null_priors(A,B,C)
%
% A{1},B{m},C  - binary constraints on extrinsic connections
%
% pE - prior expectation
% pC - prior covariance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_null_priors.m 5908 2014-03-05 20:31:57Z karl $
 
% default: a single source model
%--------------------------------------------------------------------------
if nargin < 3
    A   = {1};
    B   = {};
    C   = 1;
end

% orders
%--------------------------------------------------------------------------
n     = size(C,1);                                % number of sources
u     = size(C,2);                                % number of inputs

% parameters for Jacobian
%==========================================================================

% canonical source
%--------------------------------------------------------------------------
a     = 8;
Hz    = [8 16 48];
for i = 1:length(Hz)
   b      = 2*pi*Hz(i);
   J{i,i} = [a b;-b a];
end

% Jacobian
%--------------------------------------------------------------------------
J     = full(spm_cat(J));
nx    = length(J);
pE.A  = logm(kron(eye(n,n),J));
pC.A  = kron(A{1},ones(nx,nx));

% Bilinear terms
%--------------------------------------------------------------------------
for i = 1:u
    pE.B{i} = kron(B{i},zeros(nx,nx));
    pC.B{i} = kron(B{i},ones(nx,nx));
end

% input coeficicents
%--------------------------------------------------------------------------
pE.C  = kron(C,0);
pC.C  = C;


% input coeficicents
%--------------------------------------------------------------------------
pE.D  = kron(zeros(n,n),ones(nx,1));
pC.D  = kron(speye(n,n),ones(nx,1));

function [varargout] = spm_ind_priors(A,B,C,Nf)
% prior moments for a neural-mass model of erps
% FORMAT [pE,gE,pC,gC] = spm_ind_priors(A,B,C,dipfit,Nu,Nf)
% A{2},B{m},C  - binary constraints on extrinsic connections
% Nf           - number of frequencies

%
% pE - prior expectation - f(x,u,P,M)
% gE - prior expectation - g(x,u,G,M)
%
% spatial parameters
%--------------------------------------------------------------------------
% or gE.L    - coeficients of local modes - ECD
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A    - trial-invariant
%    pE.B{m} - trial-dependent
%    pE.C    - stimulus-stimulus dependent
%
%    p.K     - global rate [coupling] constant
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R - magnitude, onset and dispersion
%    pE.N - background fluctuations
%
% pC - prior covariances: cov(spm_vec(pE))
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

% orders
%--------------------------------------------------------------------------
n    = size(C,1);                                 % number of sources
n1   = ones(n,1);

% paramters for electromagnetic forward model
%--------------------------------------------------------------------------
G.L  = sparse(1,n);  U.L  =  sparse(1,n) + 8;

% Global scaling
%--------------------------------------------------------------------------
E.K  = 0;
V.K  = 1/4;

% set extrinsic connectivity - linear and noninlear (cross frequency)
%--------------------------------------------------------------------------
E.A  = kron(speye(Nf,Nf),A{1}/4 - speye(n,n));
V.A  = kron(speye(Nf,Nf),A{1}) + kron(1 - speye(Nf,Nf),A{2});

% input-dependent
%--------------------------------------------------------------------------
for i = 1:length(B)
    E.B{i} = sparse(n*Nf,n*Nf);
    V.B{1} = kron(ones(Nf,Nf),B{i});
end

% exognenous inputs
%--------------------------------------------------------------------------
E.C  = kron(1./[1:Nf]',C);
V.C  = E.C > 0;

% set stimulus parameters: magnitude, onset and dispersion
%--------------------------------------------------------------------------
E.R  = [0 0 1];
V.R  = [1/2 1/16 1/16];

% background fluctuations: amplitude and Hz
%--------------------------------------------------------------------------
E.N  = [0 0 10];
V.N  = [1 1 1];


% prior momments
%--------------------------------------------------------------------------
varargout{1} = E;
varargout{2} = G;
varargout{3} = diag(sparse(spm_vec(V)));
varargout{4} = diag(sparse(spm_vec(U)));
  

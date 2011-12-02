function [qP] = spm_DEM_qC(M)
% returns the conditional precision over hidden states
% FORMAT [qP] = spm_DEM_qC(M)
%
% M  - recognition  model
%   M(1).x    = Conditional expectation of hidden states
%   M(1).v    = Conditional expectation of causal states
%
% qP     - conditional precision, evaluated at M.x, M.v
%__________________________________________________________________________
%
% see spm_DEM and spm_ADEM for details.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_qC.m 4580 2011-12-02 20:22:19Z karl $
 
 
% order parameters (d = n = 1 for static models) and checks
%==========================================================================
d    = M(1).E.d + 1;                   % embedding order of q(v)
n    = M(1).E.n + 1;                   % embedding order of q(x) (n >= d)
s    = M(1).E.s;                       % smoothness - s.d. of kernel (bins)
 
% number of states and parameters - recognition model
%--------------------------------------------------------------------------
nl   = size(M,2);                      % number of levels
nv   = sum(spm_vec(M.m));              % number of v (causal states)
nx   = sum(spm_vec(M.n));              % number of x (hidden states)
ny   = M(1).l;                         % number of y (inputs)
nc   = M(end).l;                       % number of c (prior causes)
 
% precision (R) and covariance of generalised errors
%--------------------------------------------------------------------------
iV   = spm_DEM_R(n,s);
 
% precision components Q{}
%==========================================================================
Q     = {};
for i = 1:nl
    q0{i,i} = sparse(M(i).l,M(i).l);
    r0{i,i} = sparse(M(i).n,M(i).n);
end
Q0    = kron(iV,spm_cat(q0));
R0    = kron(iV,spm_cat(r0));
for i = 1:nl
    for j = 1:length(M(i).Q)
        q          = q0;
        q{i,i}     = M(i).Q{j};
        Q{end + 1} = blkdiag(kron(iV,spm_cat(q)),R0);
    end
    for j = 1:length(M(i).R)
        q          = r0;
        q{i,i}     = M(i).R{j};
        Q{end + 1} = blkdiag(Q0,kron(iV,spm_cat(q)));
    end
end
 
% and fixed components P
%--------------------------------------------------------------------------
Q0    = kron(iV,spm_cat(spm_diag({M.V})));
R0    = kron(iV,spm_cat(spm_diag({M.W})));
Qp    = blkdiag(Q0,R0);
nh    = length(Q);
 
% fixed priors on states
%==========================================================================
xP    = spm_cat(spm_diag({M.xP}));
Px    = kron(iV(1:n,1:n),speye(nx,nx)*exp(-8) + xP);
Pv    = kron(iV(1:d,1:d),speye(nv,nv)*exp(-8));
Pu    = spm_cat(spm_diag({Px Pv}));
 
% priors on parameters (in reduced parameter space)
%==========================================================================
qp.p  = cell(nl,1);
for i = 1:(nl - 1)
    
    % eigenvector reduction: p <- pE + qp.u*qp.p
    %----------------------------------------------------------------------
    qp.u{i}   = spm_svd(M(i).pC);                    % basis for parameters
    M(i).p    = size(qp.u{i},2);                     % number of qp.p
    qp.p{i}   = sparse(M(i).p,1);                    % initial qp.p
 
end
 
% initialise cell arrays for generalised states
%==========================================================================
qu.x      = cell(n,1);
qu.v      = cell(n,1);
qu.y      = cell(n,1);
qu.u      = cell(n,1);

[qu.x{:}] = deal(sparse(nx,1));
[qu.v{:}] = deal(sparse(nv,1));
[qu.y{:}] = deal(sparse(ny,1));
[qu.u{:}] = deal(sparse(nc,1));
 
% fill in values of hidden states and causes
%--------------------------------------------------------------------------
qu.x{1}   = spm_vec({M(1:end - 1).x});
qu.v{1}   = spm_vec({M(1 + 1:end).v});
 
 
% get data precisions
%==========================================================================
 
% hyperpriors
%--------------------------------------------------------------------------
qh.h  = spm_vec({M.hE M.gE});
iS    = Qp;
for i = 1:nh
    iS = iS + Q{i}*exp(qh.h(i));
end
 
% get error gradients
%--------------------------------------------------------------------------
[E dE] = spm_DEM_eval(M,qu,qp);
 
 
% conditional covariance of states
%==========================================================================
qP     = dE.du'*iS*dE.du + Pu;

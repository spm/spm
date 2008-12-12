function model = spm_mvb_G(X,L,X0,G,V)
% Multivariate Bayesian inversion of a linear model
% FORMAT model = spm_mvb_G(X,L,X0,G,V);
% X      - contrast or target vector
% L      - pattern matrix (n x m)
% X0     - confounds
% G      - pattern subsets (in columns of G) (m x h)
% V      - cell array of observation noise covariance components
%
% returns model:
%                F: log-evidence [F(0), F(1),...]
%                G: pattern switches
%                h: covariance hyperparameters (on R and cov(E))
%               qE: conditional expectation of pattern-weights
%              MAP: MAP projector (pattern-weights)
%               Cp: prior covariance (pattern space)
%__________________________________________________________________________
%
% model: X = L*P + X0*Q + R
%        P = E;           
%   cov(E) = h1*diag(G(:,1)) + h2*diag(G(:,2)) + ...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_G.m 2559 2008-12-12 17:10:23Z karl $
 
% defaults
%--------------------------------------------------------------------------
Nx    = size(X,1);
Nk    = size(G,1);
Np    = size(G,2);
if isempty(X0), X0 = zeros(Nx,1); end
try, V; catch,  V  = speye(Nx);   end
 
% null space of confounds
%--------------------------------------------------------------------------
X0    = full(X0);
R     = spm_svd(speye(size(X0,1)) - X0*pinv(X0));
L     = R'*L;
X     = R'*X;
Nx    = size(X,1);
 
% random effects (and serial correlations)
%--------------------------------------------------------------------------
if iscell(V)
    Ne    = length(V);
    for i = 1:Ne
        Qe{i} = R'*V{i}*R;
    end
elseif isnumeric(V)
    Ne = 1;
    Qe = {R'*V*R};
else
    Ne = 1;
    Qe = {speye(R'*R)};
end
 
% assemble empirical priors
%==========================================================================
Qp    = {};
LQpL  = {};
for i = 1:Np
    j       = find(G(:,i));
    Qp{i}   = sparse(j,j,1,Nk,Nk);
    LQpL{i} = L*Qp{i}*L';
end
 
% inverse solution
%==========================================================================
 
% Covariance components (with the first error Qe{1} fixed at -4)
%--------------------------------------------------------------------------
if size(L,2)
    Q = {Qe{1} Qe{:} LQpL{:}};
else
    Q = {Qe{1} Qe{:}};
end
m     = length(Q);
hE    = -32*ones(m,1);
hC    = 256*speye(m,m);
hE(1) = -4;
hC(1) = 0;
 
% ReML
%--------------------------------------------------------------------------
[Cy,h,P,F] = spm_reml_sc(X*X',[],Q,size(X,2),hE,hC);
 
h(2)  = h(2) + h(1);
h     = h(2:end);
 
% prior covariance: source space
%--------------------------------------------------------------------------
Cp    = sparse(Nk,Nk);
hp    = h([1:Np] + Ne);
for i = 1:Np
    Cp  = Cp  + hp(i)*Qp{i};
end
 
% MAP estimates of instantaneous sources
%==========================================================================
MAP   = Cp*L'*inv(Cy)*R';
qE    = MAP*R*X;
 
 
% assemble results (pattern weights)
%==========================================================================
model.F   = F;
model.G   = G;
model.h   = h;
model.qE  = qE;
model.MAP = MAP;
model.Cp  = Cp;

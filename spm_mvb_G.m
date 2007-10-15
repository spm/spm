function model = spm_mvb_G(X,Y,X0,U,G,V);
% Multivariate Bayesian inversion of a linear model
% FORMAT model = spm_mvb_G(X,Y,X0,U,G,V);
% X      - contrast or target vector
% Y      - data feature matrix
% X0     - confounds
% U      - patterns (n x m)
% G      - pattern subsets (in columns of G) (n x m)
% V      - cell array of observation noise covariance components
%
% returns model:
%                F: log-evidence [F(0), F(1),...]
%                G: pattern switches
%                U: patterns
%               qE: conditional expectation of pattern-weights
%               qC: conditional variance of pattern-weights
%                h: covariance hyperparameters
%__________________________________________________________________________
%
% model: X = Y*P + X0*Q + R
%        P = U*E;           
%   cov(E) = h1*diag(G(:,1)) + h2*diag(G(:,2)) + ...
%__________________________________________________________________________


% defaults
%--------------------------------------------------------------------------
if isempty(X0), X0 = zeros(size(X,1),1); end
if nargin < 6,  V  = speye(size(X,1));   end

% null space of confounds
%--------------------------------------------------------------------------
X0     = full(X0);
R      = orth(speye(size(X0,1)) - X0*pinv(X0));
Y      = R'*Y;
X      = R'*X;
 
% residual forming matrix
%--------------------------------------------------------------------------
ns     = size(Y,1);
nv     = size(Y,2);
nk     = size(U,2);
 
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
Np    = size(G,2);
L     = Y*U;
Qp    = {};
LQp   = {};
LQpL  = {};
for i = 1:Np
    Qp{i}   = sparse(diag(G(:,i)));
    LQp{i}  = L*Qp{i};
    LQpL{i} = LQp{i}*L';
end
 
% Inverse solution
%==========================================================================
 
% ReML
%--------------------------------------------------------------------------
if size(U,2)
    Q  = {Qe{:} LQpL{:}};
else
    Q  = Qe;
end
[Cy,h,P,F] = spm_reml_sc(X*X',[],Q,size(X,2));
 
 
% Covariances: target space - Ce and source space - L*Cp
%--------------------------------------------------------------------------
LCp   = sparse(ns,nk);
Ce    = sparse(ns,ns);
Cp    = sparse(nk,nk);
he    = h([1:Ne]);
hp    = h([1:Np] + Ne);
for i = 1:Ne
    Ce = Ce + he(i)*Qe{i};
end
for i = 1:Np
    Cp  = Cp  + hp(i)*Qp{i};
    LCp = LCp + hp(i)*LQp{i};
end
 
% MAP estimates of instantaneous sources
%==========================================================================
MAP   = LCp'*inv(Cy);
qE    = MAP*X;
 
% conditional covariance
% qC  = Cp - Cp*L'*iC*L*Cp;
%--------------------------------------------------------------------------
if nk
    R  = speye(nk,nk)/exp(16);
    qC = inv(L'*inv(Ce)*L + inv(Cp + R));
else
    qC = Cp - MAP*LCp;
end

% assemble results (pattern wieghts)
%==========================================================================
model.F   = F;
model.G   = G;
model.U   = U;
model.h   = h;
model.qE  = qE;
model.qC  = qC;
model.MAP = MAP;

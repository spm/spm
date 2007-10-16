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
%               qE: conditional expectation of pattern-weights
%               Cp: emprical prior covarisnce on pattern-weights
%                h: covariance hyperparameters (or R and cov(E)
%              MAP: MAP projector (pattern-weights)
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
 
% orders
%--------------------------------------------------------------------------
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
Nk    = size(G,1);
Np    = size(G,2);
L     = Y*U;
Qp    = {};
LQpL  = {};
for i = 1:Np
    j       = find(G(:,i));
    Qp{i}   = sparse(j,j,1,Nk,Nk);
    LQpL{i} = L*Qp{i}*L';
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
 
 
% Covariances: source space
%--------------------------------------------------------------------------
Cp    = sparse(Nk,Nk);
hp    = h([1:Np] + Ne);
for i = 1:Np
    Cp  = Cp  + hp(i)*Qp{i};
end
 
% MAP estimates of instantaneous sources
%==========================================================================
MAP   = Cp*L'*inv(Cy);
qE    = MAP*X;
 

% assemble results (pattern wieghts)
%==========================================================================
model.F   = F;
model.G   = G;
model.h   = h;
model.qE  = qE;
model.Cp  = Cp;
model.MAP = MAP;

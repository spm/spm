function model = spm_mvb(X,Y,X0,U,V,nG,sG)
% Multivariate Bayesian inversion of a linear model
% FORMAT model = spm_mvb(X,Y,X0,U,V,nG,sG)
% X      - contrast or target vector
% Y      - date feature matrix
% X0     - confounds
% U      - patterns
% V      - observation noise covariance
% nG     - number of Greedy iterations (nG = 1 => uniform hyperpriors)
%        - if not specified, the search will terminate when F falls
% sG     - size of successive subdivisions [default is 1/2)
%
% returns model:
%                F: log-evidence [F(0), F(1),...]
%                G: covariance partition indices
%                h: covariance hyperparameters
%                U: patterns
%               qE: conditional expectation of weights
%               qC: conditional variance of weights
%               Cp: prior covariance (pattern space)
%__________________________________________________________________________
%
% model: X = Y*P + X0*Q + R
%        P = U*E;           
%   cov(E) = h1*diag(G(:,1)) + h2*diag(G(:,2)) + ...
%__________________________________________________________________________
 
% defaults
%--------------------------------------------------------------------------
try, V;          catch, V  = [];        end
try, nG; aG = 0; catch, nG = 8; aG = 1; end
try, sG;         catch, sG = 1/2;       end
 
% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                 % number of samples
nv     = size(Y,2);                 % number of voxels
nk     = size(U,2);                 % number of parameters
 
% confounds
%--------------------------------------------------------------------------
if ~length(X0); X0 = zeros(ns,1);  end
if ~length(U);  U  = zeros(nv,0);  end
if ~length(V);  V  = speye(ns,ns); end
 
 
% null model
%--------------------------------------------------------------------------
model  = spm_mvb_G(X,Y,X0,sparse(nv,0),[],V);
if ~nk, return, end
 
% Initialise G and greedy search
%==========================================================================
F  = model.F;
G  = ones(nk,1);
for  i = 1:nG
 
    % invert
    %----------------------------------------------------------------------
    M  = spm_mvb_G(X,Y,X0,U,G,V);
    
    % record conditional estimates (and terminate automatically if aG)
    %----------------------------------------------------------------------
    if M.F > max(F) || i == 1
        model    = M;
    elseif aG
        break
    end
    
    % record free energy
    %----------------------------------------------------------------------
    F(i + 1)     = M.F
    
    % create new spatial support
    %----------------------------------------------------------------------
    g            = find(G(:,end));
    ng           = ceil(length(g)*sG);
    qE           = spm_en(M.qE(g,:));
    [q j]        = sort(-sum(qE.*2,2));
    q            = g(j(1:ng));
    G(q,end + 1) = 1;
    
    % break if cluster is one
    %----------------------------------------------------------------------
    if ng < 1/(1 - sG), break, end
 
end
 
% project pattern weights to feature (voxel) weights
%==========================================================================
model.F  = F;
model.U  = U;

% remove some patterns if there are too many
%--------------------------------------------------------------------------
qE       = sum(abs(model.qE),2);
[i j]    = sort(-qE);
try
    i    = j(1:2^11);
catch
    i    = j;
end
U        = model.U(:,i);
Cp       = model.Cp(i,i);
MAP      = U*model.MAP(i,:);
qE       = U*model.qE(i,:);

% conditional covariance: qC  = Cp - Cp*L'*iC*L*Cp;
%--------------------------------------------------------------------------
X0       = full(X0);
R        = orth(speye(size(X0,1)) - X0*pinv(X0));
L        = R'*Y*U;
LCp      = L*Cp;
qC       = sum((U*Cp).*U,2) - sum((U*LCp').*MAP,2);
model.M  = MAP;
model.qE = qE;                                     % conditional expecation
model.qC = max(qC,exp(-16));                       % conditional variance



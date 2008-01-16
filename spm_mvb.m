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
% Copyright (C) 2006 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb.m 1102 2008-01-16 20:54:31Z christophe $
 
% defaults
%--------------------------------------------------------------------------
try, V;          catch, V  = [];        end
try, nG; aG = 0; catch, nG = 8; aG = 1; end
try, sG;         catch, sG = 1/2;       end
 
% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                 % number of samples
nv     = size(Y,2);                 % number of voxels
np     = size(U,2);                 % number of patterns or parameters
nh     = length(V);                 % number of error components
 
% confounds
%--------------------------------------------------------------------------
if ~length(X0); X0 = zeros(ns,1);  end
if ~length(U);  U  = zeros(nv,0);  end
if ~length(V);  V  = speye(ns,ns); end
 
 
% null model
%--------------------------------------------------------------------------
model  = spm_mvb_G(X,zeros(ns,0),X0,[],V);
if ~np, return, end
 
% Initialise G and greedy search
%==========================================================================
F  = model.F;
L  = Y*U;
G  = ones(np,1);
for  i = 1:nG
 
    % invert
    %----------------------------------------------------------------------
    M  = spm_mvb_G(X,L,X0,G,V);
    
    % record conditional estimates (and terminate automatically if aG)
    %----------------------------------------------------------------------
    if M.F > max(F) || i == 1
        model    = M;
    elseif aG
        break
    end
    
    % record free-energy
    %----------------------------------------------------------------------
    F(i + 1)     = M.F;
    
    disp('log evidence & hyperparameters:')
    disp('       '),disp(F),disp(log(M.h'))
    
    % eliminate redundant components
    %----------------------------------------------------------------------
    j            = find(M.h((nh + 1):(end - 1)) < exp(-16));
    G(:,j)       = [];

    % create new spatial support
    %----------------------------------------------------------------------
    g            = find(G(:,end));
    ng           = ceil(length(g)*sG);
    [q j]        = sort(-sum(M.qE(g,:).^2,2));
    q            = g(j(1:ng));
    G(q,end + 1) = 1;
    
    % break if cluster is one
    %----------------------------------------------------------------------
    if ng < 1/(1 - sG), break, end
 
end
 
% project pattern weights to feature (voxel) weights
%==========================================================================
 
% remove some patterns if there are too many
%--------------------------------------------------------------------------
qE       = sum(model.qE.^2,2);
[i j]    = sort(-qE);
try
    i    = j(1:2^11);
catch
    i    = j;
end
L        = L(:,i);
U        = U(:,i);
Cp       = model.Cp(i,i);
MAP      = U*model.MAP(i,:);
qE       = U*model.qE(i,:);
 
% remove confounds from L = Y*U
%--------------------------------------------------------------------------
L        = L - X0*pinv(full(X0))*L;
 
% conditional covariance in voxel space: qC  = U*( Cp - Cp*L'*iC*L*Cp )*U';
%--------------------------------------------------------------------------
UCp      = U*Cp;
qC       = sum(UCp.*U,2) - sum((UCp*L').*MAP,2);

model.F  = F;
model.U  = U;
model.M  = MAP;
model.qE = qE;                                     % conditional expectation
model.Cp = Cp;                                     % prior covariance
model.qC = max(qC,exp(-16));                       % conditional variance

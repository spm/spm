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
%                U: ordered patterns
%               qE: conditional expectation of voxel weights
%               qC: conditional variance of voxel weights
%               Cp: prior covariance (ordered  pattern space)
%               cp: prior covariance (original pattern space)
%__________________________________________________________________________
%
% model: X = Y*P + X0*Q + R
%        P = U*E;           
%   cov(E) = h1*diag(G(:,1)) + h2*diag(G(:,2)) + ...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb.m 2559 2008-12-12 17:10:23Z karl $
 
% defaults (use splits +/- one standard deviation by default)
%--------------------------------------------------------------------------
try, V;          catch, V  = [];             end
try, nG; aG = 0; catch, nG = 8; aG = 1;      end
try, sG;         catch, sG = 1/2;            end
 
% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                 % number of samples
nv     = size(Y,2);                 % number of voxels
np     = size(U,2);                 % number of patterns or parameters
 
% confounds
%--------------------------------------------------------------------------
if ~length(X0), X0 = zeros(ns,1);  end
if ~length(U),  U  = zeros(nv,0);  end
if ~length(V),  V  = speye(ns,ns); end
 
% number of error components
%--------------------------------------------------------------------------
if iscell(V)
    nh = length(V);
else
    nh = 1;
end
 
% null model
%--------------------------------------------------------------------------
model  = spm_mvb_G(X,zeros(ns,0),X0,[],V);
if ~np, return, end
 
% initialise G and greedy search
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
    lnh          = log(M.h');
    
    disp('log evidence & hyperparameters:')
    fprintf('% 8.2f',F-F(1)),fprintf('\n')
    fprintf('% 8.2f',full(lnh)),fprintf('\n\n')
    
    
    % eliminate redundant components
    %----------------------------------------------------------------------
    lnh          = lnh((nh + 1):end) - lnh(end);
    j            = find(lnh < -16);
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
cp       = model.Cp;
Cp       = cp(i,i);
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
model.qE = qE;                              % conditional expectation
model.Cp = Cp;                              % prior covariance (ordered)
model.cp = cp;                              % prior covariance (original)
model.qC = max(qC,exp(-16));                % conditional variance


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
%                U: ordered patterns
%               qE: conditional expectation of voxel weights
%               qC: conditional variance of voxel weights
%               Cp: prior covariance (ordered  pattern space)
%               cp: prior covariance (original pattern space)
%__________________________________________________________________________
%
% model: X = Y*P + X0*Q + R
%        P = U*E;           
%   cov(E) = h1*diag(G(:,1)) + h2*diag(G(:,2)) + ...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb.m 2559 2008-12-12 17:10:23Z karl $
 
% defaults (use splits +/- one standard deviation by default)
%--------------------------------------------------------------------------
try, V;          catch, V  = [];             end
try, nG; aG = 0; catch, nG = 8; aG = 1;      end
try, sG;         catch, sG = 1/2;            end
 
% get orders
%--------------------------------------------------------------------------
ns     = size(Y,1);                 % number of samples
nv     = size(Y,2);                 % number of voxels
np     = size(U,2);                 % number of patterns or parameters
 
% confounds
%--------------------------------------------------------------------------
if ~length(X0), X0 = zeros(ns,1);  end
if ~length(U),  U  = zeros(nv,0);  end
if ~length(V),  V  = speye(ns,ns); end
 
% number of error components
%--------------------------------------------------------------------------
if iscell(V)
    nh = length(V);
else
    nh = 1;
end
 
% null model
%--------------------------------------------------------------------------
model  = spm_mvb_G(X,zeros(ns,0),X0,[],V);
if ~np, return, end
 
% initialise G and greedy search
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
    lnh          = log(M.h');
    
    disp('log evidence & hyperparameters:')
    fprintf('% 8.2f',F-F(1)),fprintf('\n')
    fprintf('% 8.2f',full(lnh)),fprintf('\n\n')
    
    
    % eliminate redundant components
    %----------------------------------------------------------------------
    lnh          = lnh((nh + 1):end) - lnh(end);
    j            = find(lnh < -16);
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
cp       = model.Cp;
Cp       = cp(i,i);
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
model.qE = qE;                              % conditional expectation
model.Cp = Cp;                              % prior covariance (ordered)
model.cp = cp;                              % prior covariance (original)
model.qC = max(qC,exp(-16));                % conditional variance

function model = spm_mvb(X,Y,X0,U,V,nG,sG)
% Bayesian optimisation of a multivariate linear model with a greedy search
% FORMAT model = spm_mvb(X,Y,X0,U,V,nG,sG)
%
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
%                M: MAP projector: qE = M*X
%               qE: conditional expectation of voxel weights 
%               qC: conditional variance of voxel weights
%               Cp: empirical prior covariance (ordered  pattern space)
%               cp: empirical prior covariance (original pattern space)
%__________________________________________________________________________
%
% model: X = Y*P + X0*Q + R
%        P = U*E;           
%   cov(E) = h1*diag(G(:,1)) + h2*diag(G(:,2)) + ...
%
% This routine uses a multivariate Bayesian (MVB) scheme to decode or
% recognise brain states from neuroimages. It resolves the ill-posed
% many-to-one mapping, from voxel values or data features to a target
% variable, using a parametric empirical or hierarchical Bayesian model.
% This model is inverted using standard variational techniques, in this
% case Variational Laplace, to furnish the model evidence and the
% conditional density of the model's parameters. This allows one to compare
% different models or hypotheses about the mapping from functional or
% structural anatomy to perceptual and behavioural consequences (or their
% deficits). The aim of MVB is not to predict (because the outcomes are
% known) but to enable inference on different models of structure-function
% mappings; such as distributed and sparse representations. This allows one
% to optimise the model itself and produce predictions that outperform
% standard pattern classification approaches, like support vector machines.
% Technically, the model inversion and inference uses the same empirical
% Bayesian procedures developed for ill-posed inverse problems (e.g.,
% source reconstruction in EEG).
%
% CAUTION: MVB should not be used to establish a significant mapping
% between brain states and some classification or contrast vector. Its use
% is limited to comparison of different models under the assumption
% (hyperprior) that this mapping exists. To ensure the mapping exists, use
% CVA or related approaches.
%
% See spm_mvb_ui and:
%
% Bayesian decoding of brain images.
% Friston K, Chu C, Mourão-Miranda J, Hulme O, Rees G, Penny W, Ashburner J.
% Neuroimage. 2008 Jan 1;39(1):181-205
% 
% Multiple sparse priors for the M/EEG inverse problem.
% Friston K, Harrison L, Daunizeau J, Kiebel S, Phillips C, Trujillo-Barreto 
% N, Henson R, Flandin G, Mattout J.
% Neuroimage. 2008 Feb 1;39(3):1104-20.
% 
% Characterizing dynamic brain responses with fMRI: a multivariate approach.
% Friston KJ, Frith CD, Frackowiak RS, Turner R.
% Neuroimage. 1995 Jun;2(2):166-72.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb.m 7081 2017-05-27 19:36:09Z karl $
 
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
if isempty(X0), X0 = zeros(ns,1);  end
if isempty(U),  U  = zeros(nv,0);  end
if isempty(V),  V  = speye(ns,ns); end
 
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
    if i == 1
        save_model = 1;
    else
        save_model = M.F > max(F);
    end
    if  save_model
        model      = M;
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
    [q,j]        = sort(-sum(M.qE(g,:).^2,2));
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
clear M X Y
qE       = sum(model.qE.^2,2);
[i,j]    = sort(-qE);
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

% try to save conditional expectations (if there is enough memory)
%--------------------------------------------------------------------------
try
    qE   = U*model.qE(i,:);
catch
    qE   = [];
end

% remove confounds from L = Y*U
%--------------------------------------------------------------------------
L        = L - X0*pinv(full(X0))*L;
 
% conditional covariance in voxel space: qC  = U*( Cp - Cp*L'*iC*L*Cp )*U';
%--------------------------------------------------------------------------
UCp      = U*Cp;
qC       = sum(UCp.*U,2) - sum((UCp*L').*MAP,2);
 
model.F  = F;
model.U  = U;
model.M  = MAP;                             % MAP projector
model.qE = qE;                              % conditional expectation
model.Cp = Cp;                              % prior covariance (ordered)
model.cp = cp;                              % prior covariance (original)
model.qC = max(qC,exp(-16));                % conditional variance

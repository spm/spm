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
%
% See spm_mvb and:
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
% $Id: spm_mvb_G.m 3813 2010-04-07 19:21:49Z karl $
 
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
 
% Covariance components (with the first error Qe{1} fixed)
%--------------------------------------------------------------------------
if size(L,2) > 0
    Q = [Qe LQpL];
else
    Q = Qe;
end

% ReML (with mildly informative hyperpriors) and a lower bound on error
%--------------------------------------------------------------------------
N          = size(X,2);
YY         = X*X';
Qe         = exp(-2)*trace(YY)*Q{1}/trace(Q{1})/N;
hE         = -4;
hC         = 16;
[Cy,h,P,F] = spm_reml_sc(YY,[],Q,N,hE,hC,Qe);


% prior covariance: source space
%--------------------------------------------------------------------------
Cp    = sparse(Nk,Nk);
for i = 1:Np
    Cp = Cp + h(i + Ne)*Qp{i};
end
 
% MAP estimates of instantaneous sources
%==========================================================================
MAP   = Cp*L'/Cy*R';
qE    = MAP*R*X;
 
 
% assemble results (pattern weights)
%==========================================================================
model.F   = F;
model.G   = G;
model.h   = h;
model.qE  = qE;
model.MAP = MAP;
model.Cp  = Cp;


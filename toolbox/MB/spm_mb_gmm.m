function varargout = spm_mb_gmm(varargin)
%
% FORMAT varargout = spm_mb_gmm(varargin)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


if nargin>=1 && isa(varargin{1},'char')
    [varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});
else
    [varargout{1:nargout}] = GMMimages(varargin{:});
end
%==========================================================================

%==========================================================================
function [cluster,lb,Alpha] = GMMimages(X,E, cluster, priors, lkp, mu, label, varargin)
%__________________________________________________________________________
% Fit a [Bayesian] Gaussian mixture model to observed data.

lb     = struct('X',-Inf, 'P',-Inf, 'mu',-Inf, 'A',-Inf, 'Alpha',-Inf);

m      = cluster{1};
b      = cluster{2};
W      = cluster{3};
nu     = cluster{4};
gam    = cluster{5};
Alpha  = cluster{6};

m0     = priors{1};
b0     = priors{2};
W0     = priors{3};
nu0    = priors{4};
%gam0   = priors{5}; % Unused - so far
if numel(priors)>=6
    Alpha0 = priors{6};
end

if nargin<  7, label        = []; end
if isempty(label)
    Alpha  = [];
    Alpha0 = [];
end

skip         = uint64([1 1 1]);
iter_max     = 1024;
tolerance    = 1e-4;
subiter_max  = 1024;
subtolerance = 1e-4;

if nargin>= 8, skip         = varargin{1}; end
if nargin>= 9, iter_max     = varargin{2}; end
if nargin>=10, tolerance    = varargin{3}; end
if nargin>=11, subiter_max  = varargin{4}; end
if nargin>=12, subtolerance = varargin{5}; end

[M,~] = size(m);
A     = bsxfun(@times, W, reshape(nu, 1, 1, []));

% -------------------------------------------------------------------------
% Indicate observerved data in each element of the cell array
O = false(2^M-1,M);
for i=1:(2^M-1)
    O(i,:) = fliplr(logical(dec2bin(i,M)-'0')); % Careful here
end
% -------------------------------------------------------------------------
% EM loop
for em=1:iter_max

    olb = lb;

    % ---------------------------------------------------------------------
    % sub-EM algorithm to update Mean/Precision with missing data
    % . Responsibilities (E[z]) are kept fixed
    % . Missing values (E[z*h], E[z*hh']) are updated
    % . Cluster parameters (mu,b,A/W,n) are updated
    [SS0m,SS1m,SS2m,LX,H] = suffstat_missing(m,b,W,nu,gam,Alpha, lkp,mu, X,E, label, uint64(skip));

    if ~isempty(Alpha)
        H(1,:) = 0; % Assume 0 in labels means unknown
        Alpha  = H + Alpha0;
    end

    % -----------------------------------------------------------------
    % Initialise objective function
    [lb.mu,lb.A] = kl_gausswishart(m,b,W,nu, m0,b0,W0,nu0);
    lb.Alpha     = kl_dirichlet(Alpha,Alpha0);
    lb.X         = marginalsum(SS0m, SS1m, SS2m, m,b,W,nu, O);
    lb.P         = LX - lb.X;
    LB           = lb.mu + lb.A + lb.X + lb.P + lb.Alpha;

    for i=1:subiter_max

        if numel(SS0m)~=1
            subsubiter_max = 1;
        else
            subsubiter_max = 4;
        end
        for ii=1:subsubiter_max
            % -------------------------------------------------------------
            % Infer missing suffstat
            % sum{E[z]}, sum{E[z*x]}, sum{E[z*xx']}
            [SS0,SS1,SS2] = suffstat_infer(SS0m, SS1m, SS2m, m,A, O);
            % -------------------------------------------------------------
            % Update GMM
            [m,A,b,W,nu] = updateclusters(SS0, SS1, SS2, {m0,b0,W0,nu0});
        end

        [lb.mu,lb.A] = kl_gausswishart(m,b,W,nu, m0,b0,W0,nu0);
        lb.Alpha     = kl_dirichlet(Alpha,Alpha0);
        lb.X         = marginalsum(SS0m, SS1m, SS2m, m,b,W,nu, O);
        LBo          = LB;
        LB           = lb.mu + lb.A + lb.X + lb.P + lb.Alpha;
        if numel(SS0m)==1 || LB-LBo < subtolerance
            break
        end
    end

    % ---------------------------------------------------------------------
    % Update mixing proportions
    for k1=unique(lkp)
        ind      = lkp==k1;
        gam(ind) = (SS0(ind)+eps)/(sum(SS0(ind))+sum(ind)*eps);
    end

    % ---------------------------------------------------------------------
    % Check convergence
    osum = olb.P+olb.X+olb.mu+olb.A+olb.Alpha;
    nsum =  lb.P+ lb.X+ lb.mu+ lb.A+ lb.Alpha;
    if (nsum-osum) < tolerance
        break;
    end
end
cluster = {m, b, W, nu, gam, Alpha};

%%%%%%%%%%%%%%%%
%[SS0m,SS1m,SS2m,LX,H] = suffstat_missing(m,b,W,nu,gam, lkp,mu, X,E, uint64(skip),label,Alpha);
%lb.X = marginalsum(SS0m, SS1m, SS2m, m,b,W,nu, O);
%lb.P = LX - lb.X;
%LB   = lb.P+ lb.X+ lb.mu+ lb.A;
%fprintf('%g %g %g %g    %g\n', lb.mu, lb.A, lb.X, lb.P, LB);
%%%%%%%%%%%%%%%%
%==========================================================================

%==========================================================================
function [SS0,SS1,SS2,ll,H] = suffstat_missing(m,b,W,nu,gam,Alpha, lkp,mu, X,E, label, skip)
% FORMAT [{SS0},{SS1},{SS2},ll] = suffstat_missing(m,b,W,nu,gam, lkp,mu, X,E, skip)
% Compute sufficient statistics for each missing pattern.

lnP = bsxfun(@minus,psi(Alpha),psi(sum(Alpha,1)));
if ~isempty(lnP), lnP(1,:) = 0; end
[s0b,s1b,s2b,ll,H] = spm_gmmlib('moments',double(m),double(b),double(W),double(nu),double(gam),...
                                    uint64(lkp),mu, X,E, uint64(skip), label, lnP);
H(1,:) = 0;

% Sort out the results into cell arrays
[P,K] = size(m);
SS0 = cell(1,2^P-1);
SS1 = cell(1,2^P-1);
SS2 = cell(1,2^P-1);
i0  = K; % skip the 0th element
i1  = 0;
i2  = 0;
for i=1:(2^P-1)
    o         = logical(dec2bin(i,P)-'0');
    Po        = sum(o);
    SS2{i}    = zeros(Po,Po,K);
    SS1{i}    = zeros(Po,K);
    SS0{i}    = zeros(1,K);
    SS0{i}(:) = s0b(i0+(1:K));
    SS1{i}(:) = s1b(i1+(1:(Po*K)));
    SS2{i}(:) = s2b(i2+(1:(Po*Po*K)));
    i0       = i0 + K;
    i1       = i1 + Po*K;
    i2       = i2 + Po*Po*K;
end
%==========================================================================

%==========================================================================
function [SS0,SS1,SS2] = suffstat_infer(lSS0, lSS1, lSS2, mu,A, O)
% FORMAT [SS0,SS1,SS2] = suffstat_infer(SS0, SS1, SS2, mu,A, O)
%
% SS0  - {1xK}     Oth order sufficient statistics
% SS1  - {PoxK}    1st order sufficient statistics
% SS2  - {PoxPoxK} 2nd order sufficient statistics
% mu   - PxK       Clusters' mean
% A    - PxPxK     Clusters' precision matrix
% L    - MxP       Mask of observed channels
%
% Compute expected 1st/2nd order statistics.


%--------------------------------------------------------------------------
% Dimensions
[P,K] = size(mu);

SS0 = zeros(1,K,  'double');
SS1 = zeros(P,K,  'double');
SS2 = zeros(P,P,K,'double');

for k=1:K
    ss1 = SS1(:,k);
    Ak  = A(:,:,k);
    ss2 = SS2(:,:,k);

    for i=1:size(O,1)
        io = O(i,:);                % Observed channels
        Po = sum(io);               % Number of observed channels
        if Po == 0, continue; end
        im = ~io;                   % Missing channels

        % -----------------------------------------------------------------
        % 0th order moment
        SS0k    = lSS0{i}(k);
        SS0(k)  = SS0(k) + SS0k;

        % -----------------------------------------------------------------
        % 1st order moment
        SS1k    = lSS1{i}(:,k);

        % 0) precompute stuff
        ss1o    = ss1(io,:);
        ss1m    = ss1(im,:);
        muo     = mu(io,k);
        mum     = mu(im,k);

        iAkmm   = Ak(im,im);
        iAkmm   = iAkmm + eye(size(iAkmm))*(1e-6*max(diag(iAkmm))+1e-30);
        iAkmm   = inv(iAkmm);
        SA      = iAkmm*Ak(im,io);

        % 1) observed
        ss1o    = ss1o + SS1k;

        % 2) missing
        % > t = mu(m) + A(m,m) \ A(m,o) * (mu(o) - g)
        ss1m    = ss1m + SS0k * mum;
        ss1m    = ss1m + SA * (SS0k * muo - SS1k);
        ss1(io) = ss1o;
        ss1(im) = ss1m;

        % -----------------------------------------------------------------
        % 2nd order moment: quadratic terms
        SS2k = lSS2{i}(:,:,k);

        % 0) precompute stuff
        ss2oo = ss2(io,io);
        ss2mo = ss2(im,io);
        ss2mm = ss2(im,im);
        mumum = SS0k * (mum * mum.');
        mumuo = SS0k * (muo * muo.');
        Gmuo  = SS1k * muo.';
        Gmum  = SS1k * mum.';

        % 1) observed x observed
        ss2oo = ss2oo + SS2k;

        % 2) missing x observed
        ss2mo = ss2mo + Gmum.' + SA * (Gmuo.' - SS2k);

        % 3) missing x missing
        tmp   = SA * (SS0k * muo - SS1k) * mum.';
        ss2mm = ss2mm + mumum ...
                      + (tmp+tmp') ...
                      + SA * (SS2k + mumuo - Gmuo.' - Gmuo) * SA.';

        % 4) uncertainty ~ missing
        ss2mm = ss2mm + SS0k*iAkmm;

        ss2(io,io) = ss2oo;
        ss2(im,io) = ss2mo;
        ss2(io,im) = ss2mo.';
        ss2(im,im) = ss2mm;
    end

    SS1(:,k)   = ss1;
    SS2(:,:,k) = ss2;

end
%==========================================================================

%==========================================================================
function [m,A,b,W,nu] = updateclusters(SS0,SS1,SS2,pr)
% FORMAT [m,A,b,W,nu] = updateclusters(SS0,SS1,SS2,{m0,b0,W0,nu0})
% SS0 - 0th order sufficient statistics (sum Z_i)
% SS1 - 1st order sufficient statistics (sum Z_i * X_i)
% SS2 - 2nd order sufficient statistics (sum Z_i * (X_i * X_i'))
% pr  - List of prior Gauss-Wishart parameters.
%
% Compute posterior GMM parameters from suff stats.
K   = numel(SS0);
m0  = double(pr{1});
b0  = double(pr{2});
W0  = double(pr{3});
nu0 = double(pr{4});

% -------------------------------------------------------------------------
% Mean
b  = b0  + SS0;
m  = bsxfun(@rdivide, SS1 + bsxfun(@times,b0,m0), b);

% -------------------------------------------------------------------------
% Scale/Precision
W  = zeros(size(SS2));
nu = nu0 + SS0;
for k=1:K
    W(:,:,k) = inv(SS2(:,:,k) + b0(k) * m0(:,k) * m0(:,k).' ...
                              -  b(k) *  m(:,k) *  m(:,k).' ...
                              + inv(W0(:,:,k)));
end
A  = bsxfun(@times, W, reshape(nu, [1 1 K]));
%==========================================================================

%==========================================================================
function lb = marginalsum(SS0, SS1, SS2, m,b,W,nu, O)
% [lb,const] = spm_gmm_lib('marginalsum', SS0, SS1, SS2, m,b, W,nu, O)
%
% SS0       - {1xK}   Zero-th order moment (per config)
% SS1       - {PxK}   First   order moment (per config)
% SS1       - {PxPxK} Second  order moment (per config)
% m         - PxK     Means
% W         - PxPxK   Precision/Scale matrices
% b         - 1xK     Mean degrees of freedom
% nu        - 1xK     Precision degrees of freedom
% O         - Mx1     List of existing codes
%
% lb        -         Sum of (expected) marginal likelihoods
%
% Compute the expected log-likelihood of each observation belonging to each
% cluster: lb = sum_{i,k} E[z_ik] E[ln p(g(i) | mu_k,A_k)]

P  = size(m,1);
K  = size(m,2);
lb = 0;

for i=1:size(O,1) % For each combination of missing voxels

    % ---------------------------------------------------------------------
    % Get mask of missing values and modalities (with this particular code)
    io = O(i,:);                % Observed channels
    Po = sum(io);               % Number of observed channels
    if Po == 0, continue; end
    im = ~io;                   % Missing channels
    Pm = P-Po;                  % Number of missing channels

    for k=1:K
        % /!\ Sub-covariance is different from the inverse sub-precision
        %   inv(S(o,o)) ~ W(W(o,o) - W(o,m)*W(m,m)\W(m,o), nu - Pm)
        nuo = nu(k) - Pm;
        Wo  = W(io,io,k) - W(io,im,k)*(W(im,im,k)\W(im,io,k));
        mo  = m(io,k);
        c   = - 0.5 * Po * log(2*pi) ...
              + 0.5 * wishart_elogdet(Wo,nuo) ...
              - 0.5 * nuo * mo.' * Wo * mo ...
              - 0.5 * Po / b(k);

        lb = lb + SS0{i}(k)*c + nuo*(SS1{i}(:,k).'*Wo*mo -  0.5*sum(sum(Wo.*SS2{i}(:,:,k))));
    end
end
%==========================================================================

%==========================================================================
function [klmu,klA] = kl_gausswishart(m,b,W,nu,m0,b0,W0,nu0)
[P,K] = size(m);
A     = bsxfun(@times, W, reshape(nu, [1 1 K]));

% Lower bound
klmu = 0;
klA  = 0;
for k=1:K
    ld   = wishart_elogdet(W(:,:,k),nu(k));
    % prior
    klmu = klmu - P*log(2*pi) + P*log(b0(k)) + ld ...
                - b0(k)*(m(:,k)-m0(:,k)).'*A(:,:,k)*(m(:,k)-m0(:,k)) ...
                - P*b0(k)/b(k);
    % posterior
    klmu = klmu + P*log(2*pi) - P*log(b(k)) - ld ...
                + P;
    klA  = klA - wishart_kl(W(:,:,k), nu(k), W0(:,:,k), nu0(k));
end
klmu = 0.5 * klmu;
%==========================================================================

%==== logdet ==============================================================
function ld = logdet(A)
% A  - A postive-definite square matrix
% ld - Logarithm of determinant of A
%
% Log-determinant of a positive-definite matrix.
% Cholesky factorisation is used to compute a more stable log-determinant.


% Cholseki decomposition of A (A = C' * C, with C upper-triangular)
[C, p] = chol(A);

if p > 0
    % A should usually be positive definite, but check anyway.
    warning(['Attempting to compute log determinant of matrix ' ...
             'that is not positive definite (p=%d).'], p);
    e  = eig(A);
    ld = sum(log(e(e>0)));
else
    % Because C is triangular, |C| = prod(diag(C))
    % Hence: log|C| = sum(log(diag(C)))
    % And:   log|A| = log|C'*C| = log(|C|^2) = 2 * sum(log(diag(C)))
    ld = 2 * sum(log(diag(C)));
end
%==========================================================================

% === wishart_elogdet =====================================================
function ld = wishart_elogdet(W, nu)
K  = size(W, 1);
ld = DiGamma(0.5*nu, K) + K*log(2) + logdet(W);
%==========================================================================

% === wishart_kl ==========================================================
function kl = wishart_kl(W1,nu1, W0,nu0)
% FORMAT kl = wishart_kl(W1,nu1, W0,nu0)
K  = size(W1, 1);
kl =   0.5*nu0*(logdet(W0) - logdet(W1)) ...
     + 0.5*nu1*(trace(W0\W1) - K) ...
     + 0.5*(nu1 - nu0)*DiGamma(0.5*nu1, K) ...
     + LogGamma(0.5*nu0, K) - LogGamma(0.5*nu1, K);
%==========================================================================

% === LogGamma ============================================================
function lg = LogGamma(a, p)
if nargin < 2, p = 1; end
lg = (p*(p-1)/4)*log(pi);
for i=1:p
    lg = lg + gammaln(a + (1-i)/2);
end
%==========================================================================

% === DiGamma =============================================================
function dg = DiGamma(a, p, k)
if nargin < 3
    k = 0;
    if nargin < 2, p = 1; end
end
dg = 0;
for i=1:p
    dg = dg + psi(k, a + (1-i)/2);
end
%==========================================================================



% =========================================================================
function [GaussPrior,extras] = updatehyperpars(cluster,GaussPrior,varargin)
% FORMAT [GaussPrior,extras] = updatehyperpars(cluster,GaussPrior,varargin)
%
% REQUIRED
% cluster    - 1xS cell array where cluster{s} = {{mu,b},{W,n}}
% GaussPrior - {mu0,b0,W0,nu0}
%
% OPTIONAL
% constrained - Optimise hierarchical prior on W [false]
%
% OUTPUT
% GaussPrior - New {mu0,b0,W0,nu0}
% extras     - Struct with lower bound information, etc.
%
% Update of VB-GMM hyper-parameters (m,b,W,n).

% Parse optional arguments
%--------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'updatehyperpars';
p.addParameter('constrained',0,@islogical);
p.addParameter('b0_priors',{1e-3,1e-3});
p.parse(varargin{:});
constrained = p.Results.constrained;
b0_priors   = p.Results.b0_priors;

% Parameters
S = numel(cluster); % Number of posteriors

m0  = GaussPrior{1};
b0  = GaussPrior{2};
W0  = GaussPrior{3};
nu0 = GaussPrior{4};

N = size(m0,1);
K = size(m0,2);

% pre-allocate
LogDetW0  = zeros(size(nu0));
W         = zeros(size(W0));
p         = zeros(size(nu0));
p0        = 0;

% -------------------------------------------------------------------------
%   Gauss-Wishart "mean" parameters
% -------------------------------------------------------------------------

for k=1:K

    % ---------------------------------------------------------------------
    % Update m0 (mode, closed-form)
    Lambda   = 0;
    LambdaMu = 0;
    for s=1:S
        [m,~,W,nu] = get_posteriors(cluster,s);
        Lambda     = Lambda   + nu(k)*W(:,:,k);
        LambdaMu   = LambdaMu + nu(k)*W(:,:,k)*m(:,k);
    end
    m0(:,k) = Lambda \ LambdaMu;

    % ---------------------------------------------------------------------


    % ---------------------------------------------------------------------
    % Update b0 (mode, closed-form)
    % b0 ~ \gamma(alpha,beta)
    alpha0 = b0_priors{1};
    beta0  = b0_priors{2};
    b0(k)  = 0;
    alph   = 0.5*N*S+1;
    bet    = 0;
    for s=1:S
        [m,b,W,nu] = get_posteriors(cluster,s);
        m1         = m(:,k) - m0(:,k);
        bet        = bet + 0.5*(m1.' * (nu(k)*W(:,:,k)) * m1 + N/b(k));
    end
    b0(k) = (alph+alpha0-1)/(bet+beta0);
    % ---------------------------------------------------------------------

end


% =========================================================================
% NOT CONSTRAINED
if ~constrained

    % ---------------------------------------------------------------------
    %   Gauss-Wishart "precision" parameters
    % ---------------------------------------------------------------------

    for k=1:K

        % ---
        % Set up some constants
        sumLogDet = 0;
        sumPsi    = 0;
        Wn        = 0;
        for s=1:S
            [~,~,W,nu] = get_posteriors(cluster,s);
            sumLogDet  = sumLogDet + logdet( W(:,:,k));
            sumPsi     = sumPsi    + DiGamma(nu(k)/2, N);
            Wn         = Wn        + nu(k)*W(:,:,k);
        end
        sumLogDet = sumLogDet/S;
        sumPsi    = sumPsi/S;
        Wn        = Wn/S;
       %LogDetWn  = logdet(Wn);

        % -----------------------------------------------------------------
        % Update nu0 (mode, Gauss-Newton [convex])
        E     = inf;
        for gniter=1:1000

            % -------------------------------------------------------------
            % Update W0 (mode, closed-form)
            W0(:,:,k)   = Wn/nu0(k);
            LogDetW0(k) = logdet( W0(:,:,k));
           %LogDetW0(k) = LogDetWn - N*log(nu0(k));
            % -------------------------------------------------------------

            % ---
            % Objective function
            Eprev = E;
            E = 0.5*S*nu0(k)*( LogDetW0(k) - sumLogDet - sumPsi ) ...
                + S*LogGamma(nu0(k)/2, N);

            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = 0.5*S*( LogDetW0(k) - sumLogDet - sumPsi ...
                         + DiGamma(nu0(k)/2, N) );
            H = S/4*DiGamma(nu0(k)/2, N, 1);

            % ---
            % Update
           %nu0(k) = max(nu0(k) - H\g, N-1+N*eps);
            reg   = 0.001; % Regularise so nu0 tends towards 0
            nu0(k) = max((H+reg)\(H*nu0(k)-g), N-1+N*eps);
        end
        % -----------------------------------------------------------------

    end

    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    extras.b   = b0;
    extras.m   = m0;
    extras.n   = nu0;
    extras.W   = W0;
    extras.ldW = LogDetW0;
    extras.lb  = 0;

    GaussPrior{1} = m0;
    GaussPrior{2} = b0;
    GaussPrior{3} = W0;
    GaussPrior{4} = nu0;

% =========================================================================
% CONSTRAINED
else

    lb = -inf;
    for em=1:50
        % ---
        % Starting estimate
        if p0 == 0
            p0 = 0;
            W0 = 0;
            for k=1:K
                p0 = p0 + S*nu0(k);
                for s=1:S
                    [~,~,W,nu] = get_posteriors(cluster,s);
                    W0         = W0 + inv(nu(k)*W(:,:,k));
                end
            end
            p0 = p0/K;
            W0 = W0/K;
        end

        % -----------------------------------------------------------------
        %   Gauss-Wishart "precision" parameters
        % -----------------------------------------------------------------

        for k=1:K

            % ---
            % Set up some constants
            % > compute sum E[logdet W] and sum psi(nu/2)
            logDetW  = 0;
            psiN     = 0;
            Lambda   = 0;
            for s=1:S
                [~,~,W,nu] = get_posteriors(cluster,s);
                logDetW    = logDetW  + logdet(W(:,:,k));
                psiN       = psiN     + DiGamma(nu(k)/2, N);
                Lambda     = Lambda   + nu(k)*W(:,:,k);
            end
            logDetW  = logDetW/S;
            psiN = psiN/S;


            % -------------------------------------------------------------
            % Update nu0 (mode, Gauss-Newton [convex])
            E = inf;

            % ---------------------------------------------------------
            % Update {p,W} for W0 (posterior, closed form)
            p(k)       = p0 + S*nu0(k);
           %W(:,:,k)   = inv(inv(W0) + Lambda); % NO NEED TO DO THIS MULTIPLE TIMES
            W(:,:,k)   = (W0*Lambda + eye(size(Lambda)))\W0;

            for gniter=1:100

                % Useful values
                W0(:,:,k)   = inv(wishart_e(W(:,:,k), p(k)));
                LogDetW0(k) = -wishart_elogdet(W(:,:,k), p(k));
                % ---------------------------------------------------------

                % ---
                % Objective function
                E1 = S*nu0(k)/2 * (LogDetW0(k) - logDetW - psiN) ...
                     + S*LogGamma(nu0(k)/2, N);
                E = [E E1];

                subgain = get_gain(E);
                if subgain < 1e-6
                    % Finished
                    break
                end

                % ---
                % Gradient & Hessian
                g = S/2*(LogDetW0(k) - logDetW - psiN + DiGamma(nu0(k)/2, N));
                H = S/4 * DiGamma(nu0(k)/2, N, 1);

                % ---
                % Update
                nu0(k) = max(nu0(k) - H\g, N-1+N*eps);
            end
            % ------------------------------------------------------------

        end


        % -----------------------------------------------------------------
        %   Inverse-Wishart parameters
        % -----------------------------------------------------------------

        % ---
        % Set up some constants
        % > compute sum Logdet(psi) and sum psi(m/2)
        sumlogW = 0;
        sumPsi  = 0;
        pW      = 0;
        for k=1:K
            sumlogW = sumlogW + logdet( W(:,:,k));
            sumPsi  = sumPsi  + DiGamma(p(k)/2, N);
            pW      = pW      + p(k)*W(:,:,k);
        end
        sumlogW = sumlogW/K;
        sumPsi  = sumPsi/K;
        pW      = pW/K;


        % -----------------------------------------------------------------
        % Update p0 (mode, Gauss-Newton [convex])
        E = inf;
        for gniter=1:1000

            % -------------------------------------------------------------
            % Update W0 (closed-form)
            W0 = pW/p0;
            LogDetW0 = logdet( W0);
            % -------------------------------------------------------------

            % ---
            % Objective function
            Eprev = E;
            E = p0*K/2*( N*LogDetW0 - sumlogW - sumPsi ) + K*LogGamma(p0/2, N);
            if E == Eprev
                break;
            end

            % ---
            % Gradient & Hessian
            g = K/2*( LogDetW0 - sumlogW - sumPsi + DiGamma(p0/2, N) );
            H = K/4*DiGamma(p0/2, N, 1);

            % ---
            % Update
            p0 = max(p0 - H\g, N-1+N*eps);

        end
        % -----------------------------------------------------------------


        % ---
        % Objective function
        nlb  = 0;
        for k=1:K
            nlb  = nlb - wishart_kl(W(:,:,k), p(k), W0, p0);
        end

        lb   = [lb nlb];
        gain = get_gain(lb);

        if gain < 1e-3
            % Finished
            break
        end

    end % < "EM" loop

    % ---------------------------------------------------------------------
    %   Save results
    % ---------------------------------------------------------------------
    extras.b   = b0;
    extras.m   = m0;
    extras.nu  = nu0;
    extras.W   = W0;
    extras.ldW = LogDetW0;
    extras.W   = W;
    extras.p   = p;
    extras.W0  = W0;
    extras.p0  = p0;
    extras.lb  = 0;
    for k=1:K
        extras.lb  = extras.lb - wishart_kl(W(:,:,k), p(k), W0, p0);
    end

    GaussPrior{1} = m0;
    GaussPrior{2} = b0;
    GaussPrior{3} = W0;
    GaussPrior{4} = nu0;
end
% =========================================================================

% === Dirichlet hyperparameters ===========================================
function [alpha0,lb] = dirichlet_hyperparameters(Alpha,hp)
% maximise p(Alpha | alpha0) p(alpha0 | v0, eta0)
% FORMAT [alpha0,E] = dirichlet_hyperparameters(Alpha,hp)
% Alpha - K x N array
% hp    - hyperparameter cell array of {eta0,a0}
%         Defaults are:
%             eta0 = 1e-2;
%             a0   = 100;
%         Used to compute:
%             v0   = -eta0*(psi(a0) - psi(K*a0))*ones(K,1)

% Conjugate Dirichlet (CD) priors.
% See https://en.wikipedia.org/wiki/Dirichlet_distribution#Conjugate_prior_of_the_Dirichlet_distribution
if nargin<2
    eta0 = 1e-2;
    a0   = 100;
else
    eta0 = hp{1};
    a0   = hp{2};
end
[K,N] = size(Alpha);
v0    = -eta0*(psi(a0) - psi(K*a0))*ones(K,1);
assert(all(v0>0),                     'Improper Dirichlet hyperprior: all(v0>0) not satisfied')
assert(eta0>-1,                       'Improper Dirichlet hyperprior: eta0>-1 not satisfied')
assert(eta0<0 | sum(exp(-v0/eta0))<1, 'Improper Dirichlet hyperprior: eta0<0 | sum(exp(-v0/eta0))<1 not satisfied');


%% Symbolic workings related to \ln CD(a0 | v, eta)
% lnB  = @(a) sum(gammaln(a)) - gammaln(sum(a));
% lnCD = @(a0,v,eta)(-eta*lnB(a0) - v'*a0);
% syms v1 v2 v3 nu positive
% v  = [v1 v2 v3]';
% syms la1 la2 la3 real
% la = [la1 la2 la3]';
% a0 = exp(la);
% E  = simplify(-lnCD(a0,v,eta))
% g  = (eta*(psi(a0) - psi(sum(a0))) + v).*a0
% assert(simplify(g(1) - diff(-lnCD(a0,v,nu),la1))==0)
% H  = -eta*(a0*a0').*psi(1,sum(a0)) + diag(g + eta*a0.^2.*psi(1,a0))
% assert(simplify(H(1,2) - diff(diff(-lnCD(a0,v,nu),la1),la2),1000)==0)
% assert(simplify(H(1,1) - diff(diff(-lnCD(a0,v,nu),la1),la1),1000)==0)


% Formulate everything within the CD framework
lnP   = bsxfun(@minus, psi(Alpha), psi(sum(Alpha,1))); % E[\ln P]
v     = v0 - sum(lnP,2);   % Sufficient statistic
eta   = N  + eta0;         % Number of observations

% Starting estimates (working with log(alpha0))
la0    = log(mean(Alpha,2));
alpha0 = exp(la0);
for it=1:100

    % Gradient w.r.t. \ln \alpha_0
    g   = (eta*(psi(alpha0) - psi(sum(alpha0))) + v).*alpha0;

    % Positive definite approximation of Hessian w.r.t. \ln \alpha_0
    % Using g, instead of abs(g), gives the actual Hessian.
    H   = eta*(diag(psi(1,alpha0))-psi(1,sum(alpha0))).*(alpha0*alpha0') + diag(abs(g));
    H   = H + eye(size(H))*(max(abs(diag(H)))*1e-6);

    % Newton update
    la0o   = la0;
    la0    = la0 - H\g;
    alpha0 = exp(la0);

    if norm(la0-la0o)^2 <= norm(la0)^2*1e-9, break; end
end
alpha0 = alpha0 + eps;
lb  = [(-eta0*(sum(gammaln(alpha0)) - gammaln(sum(alpha0))) - v0'*alpha0)
       (-eta *(sum(gammaln(alpha0)) - gammaln(sum(alpha0))) - v' *alpha0)];


function [alpha0,E] = dirichlet_hyperparameters_ml(Alpha)
% THIS FUNCTION IS UNUSED - WILL BE REMOVED
% Maximise p(Alpha|alpha0), where alpha0 are Dirichlet priors and Alpha
% are hyperparameters of posterior distributions.
[K,N] = size(Alpha);
lnP   = bsxfun(@minus, psi(Alpha), psi(sum(Alpha,1))); % \ln \tilde{\pi}
ss    = sum(lnP,2);        % Sufficient statistic
la0   = zeros(K,1)-log(K); % Starting estimates (working with log(a0))
for it=1:100
    alpha0  = exp(la0);
    g   = alpha0.*(N*(psi(alpha0) - psi(sum(alpha0))) - ss);              % Gradient
    H   = N*(alpha0*alpha0').*(diag(psi(1,alpha0)) - psi(1,sum(alpha0))); % E[Hessian]
    H   = H + diag(eps+max(g,0));                                         % Actual Hessian (where g>0)
    H   = H + eye(size(H))*(max(abs(diag(H)))*1e-6);
    lao = la0;
    la0 = la0 - (H\g); % Newton update
    if norm(la0-lao)^2 <= norm(la0)^2*eps, break; end
end
alpha0 = exp(la0)+eps; % Add eps because psi(0) is minus infinity
E      = (alpha0-1)'*ss + N*(gammaln(sum(alpha0)) - sum(gammaln(alpha0))); % ln p(Alpha|alpha0);
% =========================================================================

% =========================================================================
function lb = kl_dirichlet(Alpha,Alpha0)
% Compute E_{P~Dir(Alpha)}[ln p(P|Alpha0) - ln p(P|Alpha)]
if ~isempty(Alpha) && ~isempty(Alpha0)
    % Ignore zeros
    Alpha  = Alpha(2:end,:);
    Alpha0 = Alpha0(2:end,:);
else
    lb = 0;
    return
end

% E[ln P], where P ~ Dir(Alpha)
lnP = bsxfun(@minus, psi(Alpha), psi(sum(Alpha,1)));

% Expand so older versions of MATLAB can handle it
if size(Alpha0,2)==1 && size(Alpha,2)>1
    Alpha0 = repmat(Alpha0,1,size(Alpha,2));
end

% E[ln p(P|Alpha0)] - E[ln p(P|Alpha)]
% Numerically innacurate version
%lb  = sum(gammaln(sum(Alpha0,1)) - sum(gammaln(Alpha0),1) + sum((Alpha0- 1).*lnP,1)) - ...
%      sum(gammaln(sum(Alpha ,1)) - sum(gammaln(Alpha ),1) + sum((Alpha - 1).*lnP,1));

% Re-order for better numerical accuracy
lb = sum(sum((Alpha0-Alpha).*lnP)) + sum(gammaln(sum(Alpha0 ,1)) - gammaln(sum(Alpha ,1))) ...
                                   - sum(sum(gammaln(Alpha0),1)  - sum(gammaln(Alpha),1));
% =========================================================================

% === Gain ================================================================
function gain = get_gain(vals)
% FORMAT gain = get_gain(vals)
%
% vals - A vector of values
%
% gain - Computed gain
%
% Compute gain --- usually used to determine a stopping criteria when
% optimising
%__________________________________________________________________________
vals = vals(:);
gain = abs((vals(end - 1) - vals(end)) / ...
      (max(vals(isfinite(vals))) - min(vals(isfinite(vals)))));
%==========================================================================

%==========================================================================
function [m,b,W,nu] = get_posteriors(cluster,s)
m  = cluster{s}{1};
b  = cluster{s}{2};
W  = cluster{s}{3};
nu = cluster{s}{4};
%==========================================================================


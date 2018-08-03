function [varargout] = spm_GPclass(varargin)
% Gaussian process classification
% [p,F,K,theta,bic] = spm_GPclass(XX,t,lab,cov_fun,fun_args)
% Inputs:
%     XX       - cell array of dot product matrices
%                for training and testing data
%     t        - target values for training data
%     lab      - binary array indicating which are training data
%     cov_fun  - function for building covariance matrix
%     fun_args - additional arguments for covariance function
% Outputs
%     p     - Belonging probabilities
%     F     - Log-likelihood
%     K     - Covariance matrix
%     bic   - Adjustment to log-likelihood to account for hyper-parameter estimation
%
% See Chapter 3 of:
% C. E. Rasmussen & C. K. I. Williams, Gaussian Processes for Machine Learning, the MIT Press, 2006,
% ISBN 026218253X. c 2006 Massachusetts Institute of Technology. www.GaussianProcess.org/gpml
% or Bishop (2006) "Pattern Recognition and Machine Learning"
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_GPclass.m 7387 2018-08-03 15:13:57Z john $

if nargin>1 && isa(varargin{1},'char')
    switch lower(varargin{1})
    case {'likelihood'}
        if nargin >=2, theta    = varargin{2}; end
        if nargin >=3, XX       = varargin{3}; end
        if nargin >=4, t        = varargin{4}; end
        if nargin >=5, lab      = varargin{5}; end
        if nargin >=6, cov_fun  = varargin{6}; end
        if nargin >=7, fun_args = varargin{7}; else
                       fun_args = struct;      end
        [X,theta0] = cov_fun('init',XX,fun_args,lab);
        E          = objfun(theta,t,X,cov_fun);
        ll         = -E;
        varargout  = {ll};
    case {'hessian'}
        if nargin >=2, theta    = varargin{2}; end
        if nargin >=3, XX       = varargin{3}; end
        if nargin >=4, t        = varargin{4}; end
        if nargin >=5, lab      = varargin{5}; end
        if nargin >=6, cov_fun  = varargin{6}; end
        if nargin >=7, fun_args = varargin{7}; else
                       fun_args = struct;      end
        [X,theta0] = cov_fun('init',XX,fun_args,lab);
        [H,E0]     = Hessian(t,X,cov_fun,theta);
        varargout  = {H, -E0};
    case {'bic'}
        if nargin >=2, theta    = varargin{2}; end
        if nargin >=3, XX       = varargin{3}; end
        if nargin >=4, t        = varargin{4}; end
        if nargin >=5, lab      = varargin{5}; end
        if nargin >=6, cov_fun  = varargin{6}; end
        if nargin >=7, fun_args = varargin{7}; else
                       fun_args = struct;      end
        [X,theta0] = cov_fun('init',XX,fun_args,lab);
        bic        = LaplaceAdjustment(t,X,cov_fun,theta);
        varargout  = {bic};
    end
    return
end

if nargin >=1, XX       = varargin{1}; end
if nargin >=2, t        = varargin{2}; end
if nargin >=3, lab      = varargin{3}; end
if nargin >=4, cov_fun  = varargin{4}; end
if nargin >=5, fun_args = varargin{5}; else
               fun_args = struct; end

if nargin==3
    varargout{1} = gp_pred_ep_binclass(XX,t,lab);
    return
end

% Optimise hyperparameters theta, and compute marginal log-likelihood F
[X,theta] = cov_fun('init',XX,fun_args,lab);
[theta,F] = GPtrain(t,X,cov_fun,theta);

if nargout>=5, bic = LaplaceAdjustment(t,X,cov_fun,theta); end

% Predictions of unknown targets
X         = cov_fun('init',XX,fun_args);
K         = cov_fun(theta,X);
if size(t,2)>1
    p = gp_pred_lap_multiclass(K,t,lab);
else
%   p = gp_pred_lap_binclass(K,t,lab);
    p = gp_pred_ep_binclass(K,t,lab);
end
varargout = {p, F, K, theta};
if nargout>=5, varargout{5} = bic; end

%__________________________________________________________________________
%__________________________________________________________________________
function bic = LaplaceAdjustment(t,X,cov_fun,theta)
% Compute adjustment to marginal log-likelihood to account for uncertainty
% of estimating theta
H     = Hessian(t,X,cov_fun,theta);
[V,D] = eig(H);
d     = diag(D);
lam   = 1e-6;  % Cross-check against later definition
d     = max(d,lam);
bic   = 0.5*(numel(d)*log(lam) - sum(log(d)));
%__________________________________________________________________________
%__________________________________________________________________________
function [H,E0] = Hessian(t,X,cov_fun,theta)
dtheta = 0.05;
H      = zeros(numel(theta));
E0     = objfun(theta,t,X,cov_fun);
for i=1:numel(theta)
    thetap    = theta;
    thetap(i) = thetap(i)+dtheta;
    thetan    = theta;
    thetan(i) = thetan(i)-dtheta;
    Ep        = objfun(thetap,t,X,cov_fun);
    En        = objfun(thetan,t,X,cov_fun);
    H(i,i)    = (Ep+En - 2*E0)/dtheta^2;
    for j=(i+1):numel(theta)
        thetapp    = thetap;
        thetapp(j) = thetapp(j)+dtheta;
        Epp        = objfun(thetapp,t,X,cov_fun);
        thetapn    = thetap;
        thetapn(j) = thetapn(j)-dtheta;
        Epn        = objfun(thetapn,t,X,cov_fun);
        gp         = (Epp-Epn)/2/dtheta;

        thetanp    = thetan;
        thetanp(j) = thetanp(j)+dtheta;
        Enp        = objfun(thetanp,t,X,cov_fun);
        thetann    = thetan;
        thetann(j) = thetann(j)-dtheta;
        Enn        = objfun(thetann,t,X,cov_fun);
        gn         = (Enp-Enn)/2/dtheta;
        H(i,j)     = (gp-gn)/2/dtheta;
        H(j,i)     = H(i,j);
    end
end
%__________________________________________________________________________
%__________________________________________________________________________
function [theta,ll] = GPtrain(t,X,cov_fun,theta)
% Train Gaussian Process Classifier

% Run optimisation
% Currently only uses Powell's Method from Numerical Recipes
% - needs to be implemented more efficiently
m          = numel(theta);
thetai     = eye(m);
tolsc      = ones(m,1)*0.01;
[theta,ll] = spm_powell(theta,thetai,tolsc,@objfun,t,X,cov_fun);
ll         = -ll;
%__________________________________________________________________________
%__________________________________________________________________________
function E = objfun(theta,t,X,cov_fun)
% Objective function to minimise
K     = cov_fun(theta,X);
if size(t,2)==1
%   [a,F]   = gp_lap_binclass(K,t);
    [a,b,F] = gp_ep_binclass(K,t);
else
    [f,F]   = gp_lap_multiclass(K,t);
end
lam = 1e-6; % Cross-check against earlier definition
E   = -F + 0.5*lam*(theta(:)'*theta(:)); % - 0.5*numel(theta)*log(lam/(2*pi));
%fprintf(' %.8g',theta); fprintf(': %.8g\n', E);

%__________________________________________________________________________
%__________________________________________________________________________
function [f,F] = gp_lap_binclass(K,t,f)
% Find mode for Laplace approximation for binary classifiaction.
% Entirely derived from Rasmussen & Williams
% Algorithm 3.1 (page 46).
N = numel(t);
if nargin<3, f = zeros(N,1); end
for i=1:256
    sig = 1./(1+exp(-f));
    sig = min(max(sig,eps),1-eps);
    %figure(3); plot([t'; t'+1],[sig'; sig'],'-'); drawnow
    W   = sig.*(1-sig);
    sW  = sqrt(W);
    L   = chol(eye(N) + K.*(sW*sW'));
    b   = W.*f+(t-sig);
    of  = f;
    a   = b - diag(sW)*(L\(L'\(diag(sW)*K*b)));
    f   = K*a;
    %fprintf('\t%g', -0.5*a'*f + sum(t.*log(sig)) + sum((1-t).*log(1-sig)));
    if sum((f-of).^2)<(20*eps)^2*numel(f), break; end
end
%fprintf('\n');
if nargout>1
    F = -0.5*a'*f + sum(t.*log(sig)) + sum((1-t).*log(1-sig)) - sum(log(diag(L)));
end
%__________________________________________________________________________
%__________________________________________________________________________
function p = gp_pred_lap_binclass(K,t,o,f)
% Make predions using Laplace approximation for binary classification.
% Entirely derived from Rasmussen & Williams
% Algorithm 3.2 (page 47).
N = size(t,1);
if nargin<3
    o = false(size(K,1),1);
    o(1:size(t,1)) = true;
end

if nargin<4
    f = gp_lap_binclass(K(o,o),t);
end
sig = 1./(1+exp(-f));
W   = sig.*(1-sig);
sW  = sqrt(W);
L   = chol(eye(N) + K(o,o).*(sW*sW'));
M   = L'\diag(sW);
%os = RandStream.getDefaultStream;
os = rng;

p = zeros(sum(~o),1);
j = 0;
for i=find(~o)'
    j    = j + 1;
    mu   = K(o,i)'*(t-sig);
    v    = M*K(o,i);
    vr   = K(i,i) - v'*v;

    %s    = RandStream.create('mt19937ar','seed',0);
    %RandStream.setDefaultStream(s);
    rng(0,'twister');

    r    = randn(10000,1)*sqrt(vr)+mu;
    p(j) = mean(1./(1+exp(-r)));
    % Alternative approach (from Bishop's PRML)
    % kap  = (1+pi*vr/8)^(-1/2);  % Eq. 4.154
    % p(i) = 1./(1+exp(-kap*mu)); % Eq. 4.153
end
%RandStream.setDefaultStream(os);
rng(os);

%__________________________________________________________________________
%__________________________________________________________________________
function [f,F] = gp_lap_multiclass(K,t,f)
% Find mode for Laplace approximation for multi-class classification.
% Derived mostly from Rasmussen & Williams
% Algorithm 3.3 (page 50).
[N,C] = size(t);
if nargin<3, f = zeros(N,C); end
%if norm(K)>1e8, F=-1e10; return; end

for i=1:32
    f   = f - repmat(max(f,[],2),1,size(f,2));
    sig = exp(f)+eps;
    sig = sig./repmat(sum(sig,2),1,C);
    E   = zeros(N,N,C);
    for c1=1:C
        D         = sig(:,c1);
        sD        = sqrt(D);
        L         = chol(eye(N) + K.*(sD*sD'));
        E(:,:,c1) = diag(sD)*(L\(L'\diag(sD)));
       %z(c1)     = sum(log(diag(L)));
    end
    M = chol(sum(E,3));

    b = t-sig+sig.*f;
    for c1=1:C
        for c2=1:C
            b(:,c1) = b(:,c1) - sig(:,c1).*sig(:,c2).*f(:,c2);
        end
    end

    c   = zeros(size(t));
    for c1=1:C
        c(:,c1) = E(:,:,c1)*K*b(:,c1);
    end
    tmp = M\(M'\sum(c,2));
    a   = b-c;
    for c1=1:C
        a(:,c1) = a(:,c1) + E(:,:,c1)*tmp;
    end
    of = f;
    f  = K*a;
%   fprintf('%d -> %g %g %g\n', i,-0.5*a(:)'*f(:), t(:)'*f(:), -sum(log(sum(exp(f),2)),1));
    if sum((f(:)-of(:)).^2)<(20*eps)^2*numel(f), break; end
end
if nargout>1
    % Really not sure about sum(z) as being the determinant.
    % hlogdet = sum(z);

    R  = null(ones(1,C));
    sW = sparse([],[],[],N*(C-1),N*(C-1));
    for i=1:N
        ind         = (0:(C-2))*N+i;
        P           = sig(i,:)';
        D           = diag(P);
        sW(ind,ind) = sqrtm(R'*(D-P*P')*R);
    end
    hlogdet = sum(log(diag(chol(speye(N*(C-1))+sW*kron(eye(C-1),K)*sW))));
    F       = -0.5*a(:)'*f(:) + t(:)'*f(:) - sum(log(sum(exp(f),2)),1) - hlogdet;
    %fprintf('%g %g %g\n', -0.5*a(:)'*f(:) + t(:)'*f(:) - sum(log(sum(exp(f),2)),1), -hlogdet, F);
end
%__________________________________________________________________________
%__________________________________________________________________________
function p = gp_pred_lap_multiclass(K,t,o,f)
% Predictions for Laplace approximation to multi-class classification.
% Derived mostly from Rasmussen & Williams
% Algorithm 3.4 (page 51).
[N,C] = size(t);
if nargin<3
    o = false(size(K,1),1);
    o(1:size(t,1)) = true;
end

if nargin<4
    f = gp_lap_multiclass(K(o,o),t);
end

sig = exp(f);
sig = sig./repmat(sum(sig,2)+eps,1,C);
E   = zeros(N,N,C);
for c1=1:C   
    D         = sig(:,c1);
    sD        = sqrt(D);
    L         = chol(eye(N) + K(o,o).*(sD*sD'));
    E(:,:,c1) = diag(sD)*(L\(L'\diag(sD)));
end 
M   = chol(sum(E,3));
%os  = RandStream.getDefaultStream;
os  = rng;
p   = zeros(sum(~o),C);
j   = 0;
for i=find(~o)'
    j = j + 1;

    mu = zeros(C,1);
    S  = zeros(C,C);
    for c1=1:C
        mu(c1) = (t(:,c1)-sig(:,c1))'*K(o,i);
        b      = E(:,:,c1)*K(o,i);
        c      = (M\(M'\b));
        for c2=1:C
            S(c1,c2) = K(o,i)'*E(:,:,c2)*c;
        end
        S(c1,c1) = S(c1,c1) - b'*K(o,i) + K(i,i);
    end
    %s = RandStream.create('mt19937ar','seed',0);
    %RandStream.setDefaultStream(s);
    rng(0,'twister');
    nsamp  = 10000;
    r      = sqrtm(S)*randn(C,nsamp) + repmat(mu,1,nsamp);
    r      = exp(r);
    p(j,:) = mean(r./repmat(sum(r,1),C,1),2)';
end
%RandStream.setDefaultStream(os);
rng(os)
%__________________________________________________________________________
%__________________________________________________________________________
function [nut,taut,F] = gp_ep_binclass(K,t, nut,taut)
% Expectation Propagation for binary classification
%fprintf('norm(K)=%g\n', norm(K));
N    = size(t,1);
y    = t*2-1;
if nargin<3
    nut  = zeros(N,1);
    taut = zeros(N,1)+eps;
end
Sig  = K;

%os  = RandStream.getDefaultStream;
%s   = RandStream.create('mt19937ar','seed',0);
%RandStream.setDefaultStream(s);
os  = rng;
rng(0,'twister');

for it=1:128
    prev_nut = nut;
    for i=randperm(N)
        mui  = Sig(i,:)*nut;

        % Cavity parameters (eq 3.56)
        taum = max(1/Sig(i,i) - taut(i),eps);
        num  = mui/Sig(i,i) - nut(i);

        % Mean and variance
        mum   = num/taum;
        sigm  = 1/taum;

        % Marginal moments (eq 3.58)
        z     = y(i)*mum/sqrt(1+sigm);

       %NzZh  = Npdf(z)/Ncdf(z); would be unstable
        NzZh  = sqrt(2/pi)./erfcx(-z/sqrt(2));
        mom1  = y(i)*NzZh/sqrt(1+sigm);
        mom2  = NzZh*(z+NzZh)/(1+sigm);

        % Desired moments (eq 3.59)
        tauto   = taut(i);
        taut(i) = max(mom2/max(1-mom2/taum,eps),eps);
        nut(i)  = (mom1 + mom2*num/taum)/max(1-mom2/taum,eps);
        dtaut   = taut(i) - tauto;
        s       = Sig(:,i);
        Sig     = Sig - ((dtaut/(1+dtaut*s(i)))*s)*s'; % about 75% of the CPU time
    end

    s   = sqrt(taut);             % (eq 3.66)
    L   = chol(eye(N)+K.*(s*s')); % (eq 3.67)
    V   = L'\(repmat(s,1,N).*K);  % (eq 3.68)
    Sig = K - V'*V;               % (eq 3.68)
    if sqrt(sum((nut-prev_nut).^2)/sum(nut.^2)) < 1e-6
        break;
    end
end
%RandStream.setDefaultStream(os);
rng(os);

if nargout>2
    mu   = Sig*nut;
    sig  = diag(Sig);
    z    = y.*(nut.*sig-mu)./((sig.*taut - 1).*sqrt(1 - 1./(taut-1./sig)));
    F    = sum(logNcdf(z))...                                      % 3rd term
          +sum(log(1+taut./(1./sig-taut)))/2 -sum(log(diag(L)))... % 1st & 4th terms (eq 3.73)
          +(nut'*Sig*nut)/2 -sum((taut.*mu.^2 - 2*mu.*nut + sig.*nut.^2)./(sig.*taut - 1))/2; % 2nd & 5th
end
%__________________________________________________________________________
%__________________________________________________________________________
function l = logNcdf(z)
msk     = z<0;
l       = zeros(size(z));
l(msk)  = log(erfcx(-z( msk)/sqrt(2))) - log(2)-z(msk).^2/2;
l(~msk) = log((1+erf(z(~msk)/sqrt(2)))/2);
%__________________________________________________________________________
%__________________________________________________________________________
function p = gp_pred_ep_binclass(K,t,o,nut,taut)
N = size(t,1);
if nargin<3
    o = false(size(K,1),1);
    o(1:size(t,1)) = true;
end
if nargin<4
    [nut,taut] = gp_ep_binclass(K(o,o),t);
end

ss  = sqrt(taut);                    % (eq 3.66)
L   = chol(eye(N)+K(o,o).*(ss*ss')); % (eq 3.67)
z   = ss.*(L\(L'\(ss.*(K(o,o)*nut))));
p   = zeros(sum(~o),1);
j   = 0;
for i=find(~o)'
    j    = j + 1;
    fs   = K(o,i)'*(nut-z);
    v    = L'\(ss.*K(o,i));
    vf   = K(i,i)-v'*v;
    p(j) = 0.5*(1+erf(fs./sqrt(2*(1+vf))));
end
%__________________________________________________________________________
%__________________________________________________________________________



function varargout = spm_mb_classes(varargin)
% Get tissue classes
% FORMAT [P,dat] = spm_mb_classes(dat,mu,sett)
% dat  - Data structure for a subject
% mu   - Warped template data
% sett - Settings
% P    - Updated tissue classes
%
% FORMAT [dat,P] = spm_mb_classes('update_cat',dat,mu,sett)
% FORMAT l       = spm_mb_classes('LSE0',mu,ax)
% FORMAT l       = spm_mb_classes('LSE1',mu,ax)
% FORMAT mu      = spm_mb_classes('template_k1',mu,delta)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


if isa(varargin{1},'char')
    [varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});
else
    [varargout{1:nargout}] = get_classes(varargin{:});
end
%==========================================================================

%==========================================================================
function [P,dat] = get_classes(dat,mu,sett)

% Memory hungry. Needs to be addressed later.
if isfield(dat,'delta')
    mu = add_delta(mu,dat.delta);
end

if isfield(dat.model,'cat')
    % Categorical model
    [dat,P] = update_cat(dat,mu);
elseif isfield(dat.model,'gmm')
    % GMM model

    % Expand mu to include the background class.
    mu1 = template_k1(mu);
    if sett.gmm(dat.model.gmm.pop).nit_appear >0
        [dat,P] = spm_mb_appearance('update',dat,mu1,sett);
    else
        P       = exp(bsxfun(@minus,mu1(:,:,:,1:(size(mu1,4)-1)),LSE1(mu1,4)));
    end
else
    error('This should not happen');
end
if isfield(dat,'delta') && ~isempty(dat.delta)
    [dat.delta,tmp] = update_delta(dat.delta,mu,P,sett.del_settings,sett.accel);
    dat.E(1) = dat.E(1)+tmp;
end
%==========================================================================

%==========================================================================
function [delta,dE] = update_delta(delta,mu,P,del_settings,accel)
K = size(mu,4);
L = (eye(K)-1/(K+1))*del_settings;
H = L;
g = L*delta(:);
for k=1:size(mu,3)
    [g1,H1] = gradhess1(mu(:,:,k,:),P(:,:,k,:),delta,accel);
    g = g + double(reshape(sum(sum(g1,1),2),[K 1]));
    H = H + double(reshape(sum(sum(H1,1),2),[K K]));
end
dE       = 0.5*delta(:)'*L*delta(:);
delta(:) = delta(:) - H\g;
%==========================================================================

%==========================================================================
function [g,H] = gradhess1(mu,P,delta,accel)
dm    = size(mu);
K     = size(mu,4);
Ab    = 0.5*(eye(K)-1/(K+1)); % Bohnings bound on the Hessian
if nargin>=3 && ~isempty(delta)
    delta = reshape(delta,[1 1 1 K]);
    mu    = bsxfun(@plus,mu,delta);
end
H     = zeros([dm(1:3),K,K]);
g     = zeros([dm(1:3),K,1]);
sig   = softmax0(mu);
msk   = ~(all(isfinite(sig),4) & all(isfinite(P),4));
for k=1:K
    sig_k        = sig(:,:,:,k);
    tmp          = sig_k - P(:,:,:,k);
    tmp(msk)     = 0;
    g(:,:,:,k)   = tmp;
    tmp          = (sig_k - sig_k.^2)*accel + (1-accel)*Ab(k,k);
    tmp(msk)     = 0;
    H(:,:,:,k,k) = tmp;
    for k1=(k+1):K
        tmp           = (-sig_k.*sig(:,:,:,k1))*accel + (1-accel)*Ab(k,k1);
        tmp(msk)      = 0;
        H(:,:,:,k,k1) = tmp;
        H(:,:,:,k1,k) = tmp;
    end
end
%==========================================================================

%==========================================================================
function P = softmax0(mu,ax)
% safe softmax function (matches LSE0)

if nargin<2, ax = 4; end
mx  = max(mu,[],ax);
E   = exp(bsxfun(@minus,mu,mx));
den = sum(E,ax)+exp(-mx);
P   = bsxfun(@rdivide,E,den);
%==========================================================================

%==========================================================================
function [dat,P] = update_cat(dat,mu)
% Categorical model
P  = spm_mb_io('get_data',dat.model.cat.f);
sk = dat.samp;
P  = P(1:sk(1):end,1:sk(2):end,1:sk(3):end,:);

% Compute subject-specific categorical cross-entropy loss between
% segmentation and template
msk      = all(isfinite(P),4) & all(isfinite(mu),4);
tmp      =  sum(P.*mu,4) - LSE0(mu,4);
dat.E(1) = -sum(tmp(msk(:)));
dat.nvox =  sum(msk(:));
%==========================================================================

%==========================================================================
function mu1 = add_delta(mu,delta)
if isempty(delta)
    mu1 = mu;
else
    mu1 = bsxfun(@plus,mu,reshape(delta,[1 1 1 size(mu,4)]));
end
%==========================================================================

%==========================================================================
function l = LSE0(mu,ax)
% Strictly convex log-sum-exp function
% https://en.wikipedia.org/wiki/LogSumExp#A_strictly_convex_log-sum-exp_type_function
if nargin<2, ax = 4; end
mx = max(max(mu,[],ax),0);
l  = log(exp(-mx) + sum(exp(bsxfun(@minus,mu,mx)),ax)) + mx;
%==========================================================================

%==========================================================================
function mu1 = template_k1(mu,delta)
% Expand a template to include the implicit background class
ax   = 4;
if nargin>=2
    mu1  = add_delta(mu,delta);
else
    mu1  = mu;
end
lse  = LSE0(mu1,ax);
mu1  = cat(ax,bsxfun(@minus,mu1,lse), -lse);
%==========================================================================

%==========================================================================
function l = LSE1(mu,ax)
% log-sum-exp (including final class) function
if nargin<2, ax = 4; end
mx = max(mu,[],ax);
l  = log(sum(exp(bsxfun(@minus,mu,mx)),ax)) + mx;
%==========================================================================

%==========================================================================

function [mod,Z,V] = fil_fit(F,J,sett,ind,p,varargin)
% Bohning bound CCA stuff
% FORMAT [mod,Z,V] = fil_fit(F,sett,ind,p,mod,Z,Z0,P0)
% F{l}      - Nvox x M x N
% ind       - N x L
% p         - N x 1
% mod(l).mu - Nvox x M
% mod(l).W  - Nvox x M x K
% Z         - K x N
% Z0        - K x N
% P0        - K x K
%__________________________________________________________________________

% John Ashburner & Yu Yan
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


if nargin<3, error('Incorrect usage.'); end
if nargin<4 || isempty(ind), ind = true(size(F{1},3), numel(F)); end
N = size(ind,1);
if nargin<5 || isempty(p), p = ones([N,1],'single'); end

%%
% Figure out the various desired settings.
if ~isfield(sett,'v0'),  sett.v0  =  1; end
if ~isfield(sett,'b0'),  sett.b0  =  1; end
if ~isfield(sett,'nit'), sett.nit =  5; end
if ~isfield(sett,'K'),   sett.K   = 16; end

if nargin<6 || isempty(varargin{1})
    if isa(F,'cell')
        c = cell(1,numel(F));
    else
        c = cell(1);
        F = {F};
    end
    K   = sett.K;
    mod = struct('mu',c,'W',c);
    randn('seed',0);
    for l=1:numel(mod)
        if sum(ind(:,l),1)~=size(F{l},3)
            error(sprintf('Incompatible dimensions (%d ~= %d).', sum(ind(:,l),1), size(F{l},3)));
        end
        mod(l).mu = zeros(size(F{l},1), size(F{l},2),'single');
        mod(l).W  = randn(size(F{l},1), size(F{l},2), K,'single');
    end
    B0  = eye(K)*sett.b0;
else
    mod = varargin{1};
    K   = size(mod(1).W,3);
    B0  = eye(K)*sett.b0;
end
if nargin<7
    Z   = zeros(K,N,'single');
else
    Z   = varargin{2};
end
if nargin<8, Z0 = zeros(K,N);     else, Z0 = varargin{3}; end
if nargin<9, P0 = eye(K)/sett.v0; else, P0 = varargin{4}; end

%%
% Various assertions.
assert(size(Z,2)==N && size(Z0,2)==N && size(ind,1))
assert(size(Z,1)==K && size(Z0,1)==K)
assert(numel(mod)==size(ind,2))
for l=1:numel(mod)
    assert(size(mod(l).W,3)==K)
    assert(size(mod(l).W,1)==size(mod(l).mu,1))
    assert(size(mod(l).W,2)==size(mod(l).mu,2))
end

% Indices
csi = cumsum(ind,1);

% Assingn V to each z based on the average
Vc  = cell(numel(mod),1);
for l=1:numel(mod)
    Vc{l} = zeros([K K sum(ind(:,l))],'single');
end

%%
% Run the iterative variational Bayesian EM algorithm
for iter=1:sett.nit

    %%
    % Variational E-step
    V          = single(0);
    for n=1:N
        z  = Z(:,n);
        H  = P0;
        for l=1:numel(mod)
            if ind(n,l), H = H + HessZ(mod(l).W,J{l}(:,csi(n,l))); end
        end
        Vn = inv(H);

        g0 = P0*Z0(:,n);
        for subit=1:5
            g  = g0;
            for l=1:numel(mod)
                if ind(n,l), g = g + NumeratorZ(F{l}(:,:,csi(n,l)),J{l}(:,csi(n,l)),z,mod(l).mu,mod(l).W); end
            end
            z  = Vn*g;
        end

        Z(:,n) = z;
        V      = V + p(n)*Vn;
        for l=1:numel(mod)
            if ind(n,l)
                Vc{l}(:,:,csi(n,l)) = Vn;
            end
        end
    end

    % If there is a small amout of regularisation on the mean, centering 
    % is not recommended.
    Z  = bsxfun(@minus, Z, (Z*p)/sum(p));

    %%
    % Variational M-step
    for l=1:numel(mod)
        [mod(l).mu,mod(l).W] = UpdateW(F{l}, J{l}, Z(:,ind(:,l)), p(ind(:,l)), Vc{l}, mod(l).mu, mod(l).W, B0);
    end

end


%%
% Rotate to make ${\bf Z}{\bf Z}^T$ diagonal
if false % sett.do_orth
    ZZ      = Z*bsxfun(@times,p,Z');
    [~,~,R] = svd(ZZ); % Rotation to diagonalise ZZ
    Z       = R'*Z;    % Rotate the matrices.
   %Z0      = R'*Z0;
   %P0      = R'*P0*R;
    V       = R'*V*R;
    for l=1:numel(mod)
        Nvox     = size(F{l},1);
        M        = size(F{l},2);
        mod(l).W = reshape(reshape(mod(l).W,[Nvox*M,K])*R,[Nvox,M,K]);
    end
end

function [mod,Z,V] = Prune(mod,Z,V)
ZZ      = Z*bsxfun(@times,p,Z');
[~,~,R] = svd(ZZ); % Rotation to diagonalise ZZ
Z       = R'*Z;    % Rotate the matrices.
V       = R'*V*R;
vw      = 0;
for l=1:numel(mod)
    Nvox     = size(mod(l).W,1);
    M        = size(mod(l).W,2);
    K        = size(mod(l).W,3);
    mod(l).W = reshape(reshape(mod(l).W,[Nvox*M,K])*R,[Nvox,M,K]);
    vw       = vw + squeeze(sum(sum(mod(l).W.^2,1),2));
end
nz  = sqrt(vw.*sum(Z.^2,2)/size(Z,2));
ind = nz>0.001;
Z   = Z(ind,:);
V   = V(ind,ind);
for l=1:numel(mod)
    mod(l).W = mod(l).W(:,:,ind);
end


%% UpdateW
% Update the mean ($\bf\mu$) and basis functions ($\bf W$).
%%
% See Murphy's textbook.
%
% * Murphy K. _Machine learning: a probabilistic approach_ . Massachusetts
%   Institute of Technology. 2012:1-21.
function [mu,W] = UpdateW(F,J,Z,p,V,mu,W,B)
if isempty(mu), return; end

Nvox  = size(F,1);
M     = size(F,2);
N     = size(F,3);
K     = size(W,3);
A     = Abohning(M);
iA    = inv(A);

%%
% Update $\hat{\bf\mu}$.
nf = new_features;
for i=1:Nvox
    Vm       = iA./(J(i,:)*p+0.0001);
    Fi       = reshape(F(i,:,:),[M N]);
    ji       = J(i,:);
    Wi       = reshape(W(i,:,:),[M K]);
    mui      = mu(i,:)';
    msk      = ~isfinite(Fi);

    for subit=1:5
        if nf
            P  = SoftMax(Wi*Z + mui,1);
            R  = Fi + ji.*(A*mui - P);
        else
            P  = SoftMaxOld(bsxfun(@plus,Wi*Z, mui),1);
            R  = bsxfun(@plus,Fi, bsxfun(@times,ji,bsxfun(@minus, A*mui, P)));
        end
        R(msk) = 0;
        mui    = Vm*(R*p); % Update of mu
    end
    mu(i,:) = mui';
end

%%
% Update $\hat{\bf W}$.
for i=1:Nvox
    jp      = J(i,:)'.*p+eps; % Weights

    % E[Z*Z'] weighted with jp
    Vz0     = reshape(reshape(V,[K^2 N])*jp ,[K K]);
    if nf
        Vz1 = Z*(jp.*Z');
    else
        Vz1 = Z*bsxfun(@times,jp,Z');
    end

    % Covariance: inv(kron(Z*Z'+V,A) + kron(B,2*A)) = kron(inv(Z*Z'+V + 2*B),inv(A))
    Vw      = kron(inv(Vz1+Vz0 + 2*B),iA);

    Fi      = reshape(F(i,:,:),[M N]);
    ji      = J(i,:);
    Wi      = reshape(W(i,:,:),[M K]);
    mui     = mu(i,:)';
    msk     = ~isfinite(Fi);
    for subit=1:5
        Psi0       = Wi*Z;
        if nf
            P      = SoftMax(Psi0 + mui,1);
            R      = Fi + ji.*(A*Psi0 - P);
            R(msk) = 0;
            g      = reshape(R*(p.*Z'),[M*K,1]);
        else
            P      = SoftMaxOld(bsxfun(@plus, Psi0, mui),1);
            R      = bsxfun(@plus,Fi, bsxfun(@times,ji,bsxfun(@minus, A*Psi0, P)));
            R(msk) = 0;
            g      = reshape(R*bsxfun(@times,p,Z'),[M*K,1]);
        end
        Wi         = reshape(Vw*g,[M K]); % Update of W
    end
    W(i,:,:) = reshape(Wi,[1 M K]);
end



%% HessZ
% Compute Bohning's lower bound approximation to the Hessian used for updating
% the approximation to ${\bf z}$.
function H = HessZ(W,jac)
Nvox = size(W,1);
M    = size(W,2);
K    = size(W,3);
A    = Abohning(M);
H    = 0;
MK   = [M K];
for i=1:Nvox
    Wi = reshape(W(i,:,:),MK);
    H  = H + Wi'*(jac(i)*A)*Wi;
end

%% Abohning
% "Bohning bound": Hessian matrix replaced by a global lower bound in the Loewner ordering.
%
% ${\bf A} = \frac{1}{2}({\bf I}_M - \frac{1}{M+1})$
%%
% * Böhning D. _Multinomial logistic regression algorithm_ . Annals of the
%   institute of Statistical Mathematics. 1992 Mar 1;44(1):197-200.
function A = Abohning(M)
A  = 0.5*(eye(M)-1/(M+1));


%% NumeratorZ
% See Algorithm 21.1 of Murphy's textbook.
%%
% * Murphy K. _Machine learning: a probabilistic approach_ . Massachusetts
%   Institute of Technology. 2012:1-21.
function g = NumeratorZ(Fn,jac,z,mu,W)
if isempty(W), g = zeros([size(W,3),1],'single'); return; end
Nvox = size(Fn,1);
M    = size(Fn,2);
K    = size(W,3);
A    = Abohning(M);
Psi0 = reshape(reshape(W,[Nvox*M,K])*z,[Nvox,M]);
if new_features
    P    = SoftMax(Psi0+mu,2);
    r    = reshape(Fn+jac.*(Psi0*A-P),[1,Nvox*M]);
else
    P    = SoftMaxOld(Psi0+mu,2);
    r    = reshape(Fn+bsxfun(@times,jac,Psi0*A-P),[1,Nvox*M]);
end
g    = reshape(r*reshape(W,[Nvox*M,K]),[K,1]);


%% SoftMax
% Safe softmax over dimension $d$, which prevents over/underflow.
%
% $$p_k = \frac{\exp \psi_k}{\sum_{c=1}^K \exp \psi_c}$$
%
% With the constraint \psi_K=0
%
% $$p_k = \frac{\exp \psi_k}{1+\sum_{c=1}^{K-1} \exp \psi_c}$$
function P = SoftMax(Psi,d)
mx  = max(Psi,[],d);
E   = exp(Psi-mx);
P   = E./(sum(E,d)+exp(-mx));
%%
%%

% Version using bsxfun
function P = SoftMax0ld(Psi,d)
mx  = max(Psi,[],d);
E   = exp(bsxfun(@minus,Psi,mx));
P   = bsxfun(@rdivide, E, sum(E,d)+exp(-mx));

function nf = new_features
persistent nfp
if isempty(nfp)
    a = [1 1];
    try
        a  = a + a';
        nfp = true;
    catch
        nfp = false;
    end
end
nf = nfp;

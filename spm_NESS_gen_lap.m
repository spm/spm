function [F,S,Q,L,H,DS] = spm_NESS_gen_lap(P,M,x)
% Generate flow (f) at locations x
% FORMAT [F,S,Q,L,H,D] = spm_NESS_gen_lap(P,M,x)
% FORMAT [F,S,Q,L,H,D] = spm_NESS_gen_lap(P,M,U)
%--------------------------------------------------------------------------
% P.Qp    - polynomial coefficients for solenoidal operator
% P.Sp    - polynomial coefficients for Kernel
% P.Rp    - polynomial coefficients for mean
%
% F       - polynomial approximation to flow
% S       - negative potential (log NESS density)
% Q       - flow operator (R + G) with solenoidal and symmetric parts
% L       - correction term for derivatives of solenoidal flow
% H       - Hessian
% D       - potential gradients
%
% U = spm_ness_U(M)
%--------------------------------------------------------------------------
% M   - model specification structure
%
% U       - domain (of state space) structure
% U.x     - domain
% U.X     - sample points
% U.f     - expected flow at sample points
% U.J     - Jacobian at sample points
% U.b     - polynomial basis
% U.D     - derivative operator
% U.G     - amplitude of random fluctuations
% U.bG    - projection of flow operator (symmetric part: G)
% U.dQdp  - gradients of flow operator Q  w.r.t. flow parameters
% U.dbQdp - gradients of bQ w.r.t. flow parameters
% U.dLdp  - gradients of L w.r.t. flow parameters
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


% get basis set and gradients
%--------------------------------------------------------------------------
if ~isstruct(x)
    if ~iscell(x)
        x = num2cell(x(:)'); 
    end
    U   = spm_ness_U(M,x);
else
    U = x;
end

% dimensions and correction terms to flow operator
%==========================================================================
n    = numel(U.D);
nX   = size(U.b,1);

% sparse diagonal operator
%--------------------------------------------------------------------------
if nX > 1
    spd = @(x)sparse(1:numel(x),1:numel(x),x(:),numel(x),numel(x));
else
    spd = @(x)diag(x(:));
end

% flow operator (bQ)
%--------------------------------------------------------------------------
bQ    = 0;
for i = 1:numel(U.dbQdp)
    bQ = bQ + U.dbQdp{i}*P.Qp(i);
end

% correction term for solenoidal flow (L) and Kroneckor form of Q (Qp)
%--------------------------------------------------------------------------
Q     = cell(n,n);
L     = zeros(nX,n,'like',U.b(1));
F     = zeros(nX,n,'like',U.b(1));
for i = 1:n
    for j = 1:n
        bQij   = squeeze(bQ(i,j,:));
        Q{i,j} = spd(U.b*bQij + U.G(i,j));
        L(:,i) = L(:,i) - U.D{j}*bQij;
    end
end

% kernel for Hessian Sp
%--------------------------------------------------------------------------
for i = 1:n
    for j = i:n
        K(i,j,:) = U.b*P.Sp(:,i,j);
        K(j,i,:) = K(i,j,:);
        for k = 1:n
            DK(i,j,k,:) = U.D{k}*P.Sp(:,i,j);
            DK(j,i,k,:) = DK(i,j,k,:);
        end
    end
end

% expectation (mean) Rp
%--------------------------------------------------------------------------
E = U.b(:,1)*P.Rp;

% gradients D*S
%--------------------------------------------------------------------------
DS    = cell(n,1);
for k = 1:nX
    X          = U.X(k,:) - E(k,:);
    H(:,:,k)   = K(:,:,k)'*K(:,:,k);
    S(k,1)     = X*H(:,:,k)*X'/2;
    for j = 1:n
        DS{j}(k,1) = X*DK(:,:,j,k)'*K(:,:,k)*X' + H(j,:,k)*X';
    end
end

% predicted flow: F   = -Q*D*S - L
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        F(:,i) = F(:,i) - Q{i,j}*DS{j};
    end
    F(:,i) = F(:,i) - L(:,i);
end

return



% kernel for Hessian Sp
%--------------------------------------------------------------------------
for i = 1:n
    for j = i:n
        if i == j
            H(i,j,:) = (U.b*P.Sp(:,i,j)).^2;
            for k = 1:n
                DH(i,j,k,:) = 2*(U.D{k}*P.Sp(:,i,j)).*squeeze(H(i,j,:));
            end
        else
            H(i,j,:) = U.b*P.Sp(:,i,j);
            H(j,i,:) = H(i,j,:);
            for k = 1:n
                DH(i,j,k,:) = U.D{k}*P.Sp(:,i,j);
                DH(j,i,k,:) = DH(i,j,k,:);
            end
        end
    end
end

% gradients D*S
%--------------------------------------------------------------------------
DS    = cell(n,1);
for k = 1:nX
    X      = U.X(k,:) - E(k,:);
    S(k,1) = X*H(:,:,k)*X'/2;
    for j = 1:n
        DS{j}(k,1) = X*DH(:,:,j,k)*X' + H(j,:,k)*X';
    end
end

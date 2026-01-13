function [F,S,Q,L,H,DS,E] = spm_NESS_gen_lap(P,M,x,OPT)
% Generate flow (F) at locations x
% FORMAT [F,S,Q,L,H,D,E] = spm_NESS_gen_lap(P,M)
% FORMAT [F,S,Q,L,H,D,E] = spm_NESS_gen_lap(P,M,x)
% FORMAT [F,S,Q,L,H,D,E] = spm_NESS_gen_lap(P,M,U)
%--------------------------------------------------------------------------
% P.Qp    - polynomial coefficients for solenoidal operator
% P.Sp    - polynomial coefficients for Kernel (suprisal)
% P.Rp    - polynomial coefficients for mean   (suprisal)
%
% F       - polynomial approximation to flow
% S       - negative potential (log NESS density)
% Q       - flow operator (R + G) with solenoidal and symmetric parts
% L       - correction term for derivatives of solenoidal flow
% H       - Hessian
% D       - potential gradients
% E       - expectation of NESS density
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
% U.dLdp  - gradients of L  w.r.t. flow parameters
%
% This routine returns the flow at a number of specified points in state
% space under a polynomial approximation to any Helmholtz-Hodge
% decomposition of nonequilibrium steady-state dynamics. This routine
% shares many of the same constructs as spm_NESS_fx; with the exception
% that the dissipative part of the flow (i.e., amplitude of random
% fluctuations G) is specified as the precision of a state-space model:
% M.W.
%
% In brief, spm_NESS_fx returns the flow for a particular point in state
% space, whereas this routine returns the flow for an arbitrary number of
% points at the same time. These flows can then be used as a generative
% model for the flow sampled over grid points in state space—or during the
% evolution of some path through state space.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging

% get basis set and gradients
%--------------------------------------------------------------------------
if nargin < 3

    % assume M.X is specified
    %----------------------------------------------------------------------
    U   = spm_ness_U(M);

else

    % get basis set and gradients
    %----------------------------------------------------------------------
    if ~isstruct(x)
        if ~iscell(x)
            x = num2cell(x(:)');
        end
        U   = spm_ness_U(M,x);
    else
        U = x;
    end
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
K     = zeros(n,n,nX);
DK    = zeros(n,n,n,nX);
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
p     = 1:size(P.Rp,1);
DX    = zeros(n,n,nX);
for i = 1:n
    for j = 1:n
        if i == j
            DX(i,j,:) = 1;
        else
            %%% DX(i,j,:) = - U.D{j}(:,p)*P.Rp(p,i);
        end
    end
end

% gradients D*S
%--------------------------------------------------------------------------
DS    = cell(n,1);
DXHX  = zeros(n,1);
S     = zeros(nX,1);
H     = zeros(n,n,nX);
for k = 1:nX

    % expectation (mean) Rp and Hessian H
    %----------------------------------------------------------------------
    E        = U.b(k,p)*P.Rp;

    % gradients D*S: S  = (X - E)'*K(X)'*K(X)*(X - E)/2   =>
    %                DS =  DX'*H(X)*(X - E) + (X - E)'*DK'*K(X)*(X - E)
    %--------------------------------------------------------------------------
    X        = (U.X(k,:) - E)';
    H(:,:,k) = K(:,:,k)'*K(:,:,k);
    S(k,1)   = X'*H(:,:,k)*X/2;

    DXH   = DX(:,:,k)'*H(:,:,k);
    for i = 1:n
        DXHX(i)    = X'*DK(:,:,i,k)*K(:,:,k)*X;
        DS{i}(k,1) = DXHX(i) + DXH(i,:)*X;
    end

end

% predicted flow: F = -Q*D*S - L
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        F(:,i) = F(:,i) - Q{i,j}*DS{j};
    end
    F(:,i) = F(:,i) - L(:,i);
end

% required output
%--------------------------------------------------------------------------
if nargin > 3
    F = eval(OPT);
end



return
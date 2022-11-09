function [F,S,Q,L,H,DS] = spm_NESS_gen(P,M,U)
% Generate flow (f) at locations (U.X)
% FORMAT [F,S,Q,L,H,D] = spm_NESS_gen(P,M)
% FORMAT [F,S,Q,L,H,D] = spm_NESS_gen(P,M,U)
% FORMAT [F,S,Q,L,H,D] = spm_NESS_gen(P,M,X)
%--------------------------------------------------------------------------
% P.Qp    - polynomial coefficients for solenoidal operator
% P.Sp    - polynomial coefficients for potential
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
% Required fields:
%    M.X  - sample points
%    M.W  - (n x n) - precision matrix of random fluctuations
%    M.K  - order of polynomial expansion
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


% model specification
%--------------------------------------------------------------------------
if nargin > 2
    if ~isstruct(U), M.X = U; U = []; end
end


% use M.fs if specified
%--------------------------------------------------------------------------
if nargout < 3 && isfield(M,'fs')
    w     = {M.W(1)};
    Sp    = num2cell(P.Sp);
    Qp    = num2cell(P.Qp);
    for i = 1:size(M.X,1)
        
        % flow
        %------------------------------------------------------------------
        x      = num2cell(M.X(i,:));
        F(i,:) = M.fs(Qp{:},Sp{2:end},w{:},x{:});
        
        % negative potential
        %------------------------------------------------------------------
        if nargout == 2
            S(i)   = M.ss(Sp{:},x{:});
        end
    end
    
    return
end


% get basis or expansion from M.X (or M.x)
%--------------------------------------------------------------------------
if nargin < 3 || ~isstruct(U)
    if isfield(M,'f')
        M = rmfield(M,'f');
    end
    U     = spm_ness_U(M);
end

% dimensions and correction terms to flow operator
%==========================================================================
n  = numel(U.D);
nX = size(U.b,1);

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

% correction term for solenoidal flow (L) and Kroneckor form of Q
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


% predicted flow: F   = Q*D*S - L
%--------------------------------------------------------------------------
DS    = cell(n,1);
for j = 1:n
    DS{j} = U.D{j}*P.Sp;
end
for i = 1:n
    for j = 1:n
        F(:,i) = F(:,i) + Q{i,j}*DS{j};
    end
    F(:,i) = F(:,i) - L(:,i);
end

if nargout == 1, return, end

% (scalar) potential:  S = -log(p(x))
%--------------------------------------------------------------------------
S     = -U.b*P.Sp;

% Hessian D*D*S
%--------------------------------------------------------------------------
HH    = cell(n,n);
for i = 1:n
    for j = 1:n
       HH{i,j} = -U.H{i,j}*P.Sp;    
    end
end
H     = zeros(n,n,nX,'like',U.b(1));
for i = 1:n
    for j = 1:n
        for k = 1:nX
            H(i,j,k) = HH{i,j}(k);
        end
    end
end

function F = spm_fx_NESS(x,u,P,M)
% Generate flow (f) at locations x
% FORMAT f = spm_fx_NESS(x,u,P,M)
%--------------------------------------------------------------------------
% x       - latent states
% u       - exogenous input
% P.Qp    - polynomial coefficients for solenoidal operator
% P.Sp    - polynomial coefficients for kernel (suprisal)
% P.Rp    - polynomial coefficients for mean   (suprisal)
% P.W     - precision of random fluctuations
%
% f       - polynomial approximation to flow
% S       - negative potential (log NESS density)
% Q       - flow operator (R + G) with solenoidal and symmetric parts
% L       - correction term for derivatives of solenoidal flow
% H       - Hessian
%
% U = spm_NESS_U(x)
%--------------------------------------------------------------------------
% U.b     - polynomial basis
% U.D     - derivative operator
% U.dbQdp - gradients of bQ w.r.t. flow parameters
%
% flow  (Jacobian J = df/dx)
%--------------------------------------------------------------------------
%   S     = ln(p(x))       =  X*H*X'/2;
%   f     = (R + G)*dS/dx  =  Q*dS/dx
%   df/dx = J = -(R + G)*H = -Q*H
%   H     = K'*K
%   G     = inv(W)/2;
%
% This routine returns the flow or dynamics of a nonequilibrium
% steady-state model, parameterised under the Helmholtz-Hodge decomposition
% (up to 2nd order) in terms of flow operators Q and surprisal gradients.
% The flow has solenoidal (R) and dissipative (G) parts, while the
% surprisal S or self-information is parameterised to 2nd order, with the
% proviso that the requisite Hessian (i.e., precision of the nonequilibrium
% steady-state density) is itself state-dependent to 1st order (by
% parameterising the Hessian in terms of the outer product of a surprisal
% kernel with itself H = K’*K. The underlying surprisal or log Ness
% density, the flow operators and the Hessian can also be returned by
% increasing the number of output arguments.
%
% NB: these output arguments are suppressed in the current code because the
% accompanying inversion schemes expect the Jacobian to be returned as the
% second argument. And the Jacobian is most simply evaluated numerically,
% or with automatic differentiation.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
L = 2;
if nargin > 3
    if isfield(M,'L')
        L = M.L;
    end
end

% get (polynomial) expansion of x
%--------------------------------------------------------------------------
x   = num2cell(x(:)');
U   = spm_NESS_U(x,L);

% dimensions and correction terms to flow operator
%==========================================================================

% sparse diagonal operator
%--------------------------------------------------------------------------
n   = numel(U.D);
spd = @(x)diag(x(:));

% flow operator (bQ)
%--------------------------------------------------------------------------
bQ  = 0;
for i = 1:numel(U.dbQdp)
    bQ = bQ + U.dbQdp{i}*P.Qp(i);
end

% diffusion tensor or amplitude of random fluctuations (G)
%--------------------------------------------------------------------------
if isscalar(P.W)
    G = eye(n,n)/P.W/2;
else
    G = inv(P.W)/2;
end

% correction term for solenoidal flow (L) and Kroneckor form of Q (Qp)
%--------------------------------------------------------------------------
Q     = zeros(n,n);
L     = zeros(n,1,'like',U.b(1));
for i = 1:n
    for j = 1:n
        bQij   = squeeze(bQ(i,j,:));
        Q(i,j) = spd(U.b*bQij + G(i,j));
        L(i)   = L(i) - U.D{j}*bQij;
    end
end

% kernel for Hessian Sp
%--------------------------------------------------------------------------
K   = zeros(n,n);
DK  = zeros(n,n,n);
for i = 1:n
    for j = i:n
        K(i,j) = U.b*P.Sp(:,i,j);
        K(j,i) = K(i,j);
        for k = 1:n
            DK(i,j,k) = U.D{k}*P.Sp(:,i,j);
            DK(j,i,k) = DK(i,j,k);
        end
    end
end

% expectation (mean) Rp
%--------------------------------------------------------------------------
j     = 1:size(P.Rp,1);
E     = U.b(1,j)*P.Rp;

% gradients D*S
%--------------------------------------------------------------------------
DS    = zeros(n,1);
X     = U.X - E;
H     = K'*K;
for j = 1:n
    DS(j) = X*DK(:,:,j)'*K*X' + H(j,:)*X';
end

% predicted flow: F   = -Q*D*S - L
%--------------------------------------------------------------------------
F = -Q*DS - L;

return

function U = spm_NESS_U(x,L)
% Nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT U = spm_NESS_U(x,L)
%--------------------------------------------------------------------------
% x       - sample point
%
% U.X     - sample points
% U.b     - polynomial basis
% U.D     - derivative operator
% U.dbQdp - gradients of bQ w.r.t. flow parameters
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% state space
%--------------------------------------------------------------------------
[X,x] = spm_ndgrid(x);
[b,D] = spm_polymtx(x,L);
n     = numel(x);

% sparse diagonal operator
%--------------------------------------------------------------------------
spd = @(x)diag(x(:));

% orders of polynomial expansion
%--------------------------------------------------------------------------
nu    = (n^2 + n)/2;
nb    = size(b,2);
nB    = nu*nb;

% derivatives of flow operator Q
%--------------------------------------------------------------------------
q     = 0;
dQ    = zeros(n,n);
dbQ   = zeros(n,n,nb);
dQdp  = cell(nB,1);
dbQdp = cell(nB,1);
for i = 1:n
    for j = i:n
        for k = 1:nb
            
            % initialise partial derivatives
            %--------------------------------------------------------------
            q        = q + 1;
            dq       = dQ;
            dbQdp{q} = dbQ;
            
            if i == j
                
                % with respect to coefficients - dQdp
                %----------------------------------------------------------
                dq(i,j) = 1;
                dQdp{q} = kron(dq,spd(b(:,k)));
                
                % with respect to coefficients - dbQdp
                %----------------------------------------------------------
                dbQdp{q}(i,j,k) =  1;
                
            else
                
                % with respect to coefficients - dQdp
                %----------------------------------------------------------
                dq(i,j) =  1;
                dq(j,i) = -1;
                dQdp{q} = kron(dq,spd(b(:,k)));
                
                % with respect to coefficients - dbQdp
                %----------------------------------------------------------
                dbQdp{q}(i,j,k) =  1;
                dbQdp{q}(j,i,k) = -1;
                
            end
        end
    end
end


% assemble U structure
%--------------------------------------------------------------------------
U.X     = X;                     % sample point
U.b     = b;                     % polynomial basis
U.D     = D;                     % derivative operator
U.dbQdp = dbQdp;                 % gradients of bQ w.r.t. flow parameters



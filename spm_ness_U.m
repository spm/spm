function [U] = spm_ness_U(M,x)
% nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT [U] = spm_ness_U(M,x)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.f   - dx/dt   = f(x,u,P)  {function string or m-file}
%    M.pE  - P       = parameters of equation of motion
%    M.x   - (n x 1) = x(0) = expansion point
%    M.W   - (n x n) - precision matrix of random fluctuations
%    M.X   - sample points
%    M.K   - order of polynomial expansion
%
% U.x     - domain
% U.X     - sample points
% U.f     - expected flow at sample points
% U.J     - Jacobian at sample points
% U.b     - orthonormal polynomial basis
% U.D     - derivative operator
% U.G     - amplitude of random fluctuations
% U.v     - orthonormal operator
% U.bG    - projection of flow operator (symmetric part: G)
% U.dQdp  - gradients of flow operator Q  w.r.t. flow parameters
% U.dbQdp - gradients of bQ w.r.t. flow parameters
% U.dLdp  - gradients of L w.r.t. flow parameters
% U.K     - order polynomial expansion

%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_hd.m 8000 2020-11-03 19:04:17QDb karl $

% event space: get or create X - coordinates of evaluation grid
%--------------------------------------------------------------------------
n  = length(M.x);
if isfield(M,'FUN')
    FUN = M.FUN;
else
    FUN = 'POLY';
end

K   = 3;
if nargin < 2
    
    % use M.X
    %----------------------------------------------------------------------
    X     = M.X;
    for i = 1:size(X,1)
        x               = num2cell(X(i,:));
        [bi,Di,Hi,o]    = spm_polymtx(x,K,FUN);
        b(i,:)          = full(bi);
        for j = 1:n
            D{j,1}(i,:) = full(Di{j});
            for k = 1:n
                H{j,k}(i,:) = full(Hi{j,k});
            end
        end

    end

    % use M.X
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        x{i}  = unique(X(:,i));
    end

else
    
    % use x
    %----------------------------------------------------------------------
    [X,x]     = spm_ndgrid(x);
    [b,D,H,o] = spm_polymtx(x,K,FUN);
    
end

% constraints on flow operator
%==========================================================================
k = sum(o) < K;
c = b(:,k);
V = cell(n,1);
for i = 1:n
    V{i} = D{i}(:,k);
end


% size of subspace (nx)
%--------------------------------------------------------------------------
nx    = ones(1,n);
for i = 1:n
    nx(i) = numel(x{i});
end

% sparse diagonal operator
%--------------------------------------------------------------------------
nX      = size(X,1);
if nX > 1
    spd = @(x)sparse(1:numel(x),1:numel(x),x(:),numel(x),numel(x));
else
    spd = @(x)diag(x(:));
end

% diffusion tensor or amplitude of random fluctuations (G)
%--------------------------------------------------------------------------
if isscalar(M.W)
    G = eye(n,n)/M.W/2;
else
    G = inv(M.W)/2;
end

% flow  (Jacobian J = df/dx)
%--------------------------------------------------------------------------
%   S     = ln(p(x))
%   f     = (R + G)*dS/dx  =  Q*dS/dx
%   df/dx = J = -(R + G)*H = -Q*H
%--------------------------------------------------------------------------
X     = full(X);
f     = zeros(n,nX,'like',X);
J     = zeros(n,n,nX,'like',X);
for i = 1:nX
    
    % Jacobian at this point in state space
    %----------------------------------------------------------------------
    s        = X(i,:)';
    F        = M.f(s,0,M.pE);
    if isfield(M,'J')
        J(:,:,i) = full(M.J(s,0,M.pE,M));
    else
        J(:,:,i) = full(spm_diff(M.f,s,0,M.pE,M,1));
    end
    f(:,i)   = full(F);
end


% orthonormal polynomial expansion
%--------------------------------------------------------------------------
nu  = (n^2 - n)/2;
if size(b,1) > 1
    
    % with coefficients p: x = b*p = dxdp*p for log density Sp
    %----------------------------------------------------------------------
    v     = b\spm_orth(b,'norm'); 
    b     = b*v;
    for i = 1:n
        D{i} = D{i}*v;
        for j = 1:n
            H{i,j} = H{i,j}*v;
        end
    end
    
    % with coefficients p: x = c*p = dxdp*p for flow operators Qp
    %----------------------------------------------------------------------
    u     = c\spm_orth(c,'norm'); 
    c     = c*u;
    for i = 1:n
        V{i} = V{i}*u;
    end
    u     = kron(eye(nu,nu),u);

else
    v = 1;
    u = 1;
end
nc    = size(c,2);
nB    = nu*nc;

% coefficients (bQ) of flow operator Q (symmetric part)
%--------------------------------------------------------------------------
bG    = zeros(n,n,nc,'like',c);
for i = 1:n
    for j = 1:n
        bG(i,j,1) = G(i,j)/c(1);
    end
end

% derivatives of flow operator Q
%--------------------------------------------------------------------------
q     = 0;
dQ    = zeros(n,n);
dbQ   = zeros(n,n,nc);
dQdp  = cell(nB,1);
dbQdp = cell(nB,1);
for i = 1:n
    for j = (i + 1):n
        for k = 1:nc
            
            % with respect to coefficients - dQdp
            %--------------------------------------------------------------
            q       = q + 1;
            dq      = dQ;
            dq(i,j) =  1;
            dq(j,i) = -1;
            dQdp{q} = kron(dq,spd(c(:,k)));
            
            % with respect to coefficients - dbQdp
            %--------------------------------------------------------------
            dbQdp{q}        = dbQ;
            dbQdp{q}(i,j,k) =  1;
            dbQdp{q}(j,i,k) = -1;
        end
    end
end

% derivatives of L with respect to coefficients - dLdp
%--------------------------------------------------------------------------
dLdp  = cell(nB,1); [dLdp{:}] = deal(zeros(nX,n,'like',X));
for i = 1:n
    for j = 1:n
        for k = 1:nB
            dLdp{k}(:,i) = dLdp{k}(:,i) - V{j}*reshape(dbQdp{k}(i,j,:),nc,1);
        end
    end
end

% assemble U structure
%--------------------------------------------------------------------------
U.x     = x;                     % domain
U.X     = X;                     % sample points
U.f     = f;                     % expected flow at sample points
U.J     = J;                     % Jacobian at sample points
U.b     = b;                     % orthonormal polynomial basis (S)
U.c     = c;                     % orthonormal polynomial basis (Q)
U.D     = D;                     % derivative operator (S)
U.V     = V;                     % derivative operator (Q)
U.H     = H;                     % Hessian operator
U.G     = G;                     % amplitude of random fluctuations
U.v     = v;                     % orthonormal operator
U.u     = u;                     % orthonormal operator (Kroneckor form)
U.bG    = bG;                    % projection of flow operator
U.dQdp  = dQdp;                  % gradients of Q  w.r.t. flow parameters
U.dbQdp = dbQdp;                 % gradients of bQ w.r.t. flow parameters
U.dLdp  = dLdp;                  % gradients of L  w.r.t. flow parameters
U.nx    = nx;                    % dimensions
U.o     = o;                    % orders


return


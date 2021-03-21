function U = spm_ness_U(M,x)
% nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT U = spm_ness_U(M,x)
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
% U.nx    - dimensions
% U.o     - orders

%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_hd.m 8000 2020-11-03 19:04:17QDb karl $

% event space: get or create X - coordinates of evaluation grid
%--------------------------------------------------------------------------

if isfield(M,'FUN')
    FUN = M.FUN;
else
    FUN = 'POLY';
end
if isfield(M,'K')
    K   = M.K;
else
    K   = 3;
end

if nargin < 2
    
    % use M.X
    %----------------------------------------------------------------------
    X      = M.X;
    [nX,n] = size(X);
    for  i = 1:nX
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
    for i = 1:n
        x{i}  = unique(X(:,i));
    end

else
    
    % use x
    %----------------------------------------------------------------------
    [X,x]     = spm_ndgrid(x);
    [b,D,H,o] = spm_polymtx(x,K,FUN);
    [nX,n]    = size(X);
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
DDG = 0;
if DDG
    nu = (n^2 - n)/2;
else
    nu = (n^2 + n)/2;
end

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
    
    % Kroneckor form for flow operators Qp
    %----------------------------------------------------------------------
    u     = kron(eye(nu,nu),v);

else
    v = 1;
    u = 1;
end
nb    = size(b,2);
nB    = nu*nb;

% coefficients (bQ) of flow operator Q (symmetric part)
%--------------------------------------------------------------------------
bG    = zeros(n,n,nb,'like',b);
for i = 1:n
    for j = 1:n
        bG(i,j,1) = G(i,j)/b(1);
    end
end

% derivatives of flow operator Q
%--------------------------------------------------------------------------
q     = 0;
dQ    = zeros(n,n);
dbQ   = zeros(n,n,nb);
dQdp  = cell(nB,1);
dbQdp = cell(nB,1);
for i = 1:n
    for j = (i + DDG):n
        for k = 1:nb
            
            % initialise partial derivatives
            %----------------------------------------------------------
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

% derivatives of L with respect to coefficients - dLdp
%--------------------------------------------------------------------------
dLdp  = cell(nB,1); [dLdp{:}] = deal(zeros(nX,n,'like',X));
for i = 1:n
    for j = 1:n
        for k = 1:nB
            dLdp{k}(:,i) = dLdp{k}(:,i) - D{j}*reshape(dbQdp{k}(i,j,:),nb,1);
        end
    end
end

% assemble U structure
%--------------------------------------------------------------------------
U.x     = x;                     % domain
U.X     = X;                     % sample points
U.f     = f;                     % expected flow at sample points
U.J     = J;                     % Jacobian at sample points
U.b     = b;                     % orthonormal polynomial basis
U.D     = D;                     % derivative operator
U.H     = H;                     % Hessian operator
U.G     = G;                     % amplitude of random fluctuations
U.v     = v;                     % orthonormal operator
U.u     = u;                     % orthonormal operator (Kroneckor form)
U.bG    = bG;                    % projection of flow operator
U.dQdp  = dQdp;                  % gradients of Q  w.r.t. flow parameters
U.dbQdp = dbQdp;                 % gradients of bQ w.r.t. flow parameters
U.dLdp  = dLdp;                  % gradients of L  w.r.t. flow parameters
U.nx    = nx;                    % dimensions
U.o     = o;                     % orders


return


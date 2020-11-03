function [p0,X,F,f,H,J,E] = spm_ness_hd(M,x)
% nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT [p0,X,F,f,H,J,E] = spm_ness_hd(M,x)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.f   - dx/dt   = f(x,u,P)  {function string or m-file}
%    M.pE  - P       = parameters of equation of motion
%    M.x   - (n x 1) = x(0) = expansion point
%    M.W   - (n x n) - precision matrix of state noise
% x    - cell array of vectors specifying evaluation grid
%
% p0   - nonequilibrium steady-state
% X    - evaluation points of state space
% F    - expected flow
% f    - original flow
% H    - expected Euclidian norm of Hessian
% J    - expected Euclidian norm of Jacobian
% E    - Lyapunov exponents
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_hd.m 8000 2020-11-03 19:04:17Z karl $


% event space: get or create X - coordinates of evaluation grid
%--------------------------------------------------------------------------
n  = length(M.x);
if nargin < 2
    
    % use M.X
    %----------------------------------------------------------------------
    X     = M.X;
    for i = 1:size(X,2)
        x{i} = unique(X(:,i));
    end
else
    
    % use x
    %----------------------------------------------------------------------
    [X,x] = spm_ndgrid(x);
    
end

% size of subspace
%--------------------------------------------------------------------------
nX    = size(X,1);
nY    = nX*n;
dx    = ones(1,n);
nx    = ones(1,n);
for i = 1:n
    nx(i) = numel(x{i});
    if nx(i) > 1
        dx(i) = x{i}(2) - x{i}(1);
    end
end

% differential operators
%==========================================================================
spd = @(x)spdiags(x(:),0,numel(x),numel(x));

% gradients based upon Helmholtz decomposition
%==========================================================================

% diffusion tensor more amplitude of random fluctuations (G)
%--------------------------------------------------------------------------
if isscalar(M.W)
    G = eye(n,n)/M.W/2;
else
    G = inv(M.W)/2;
end


% flow constraints (Jacobian J = df/dx)(R = -Q')
%--------------------------------------------------------------------------
%   S     = ln(p(x))
%   f     = (R + G)*dS/dx
%   df/dx = J = -(R + G)*H
%--------------------------------------------------------------------------
f     = zeros(n,nX);
J     = zeros(n,n,nX);
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

I     = speye(n,n);
Q     = zeros(n,n,nX);
d     = speye(n*n)/512;
for i = 1:nX
    
    % Hessian and solenoidal operator
    %------------------------------------------------------------------
    Z  = kron(I,J(:,:,i)) + kron(conj(J(:,:,i)),I);
    Y  = spm_vec(J(:,:,i)*G - G*J(:,:,i)');
    R  = reshape((Z'*Z + d)\(Z'*Y),n,n);

    % potential gradients
    %------------------------------------------------------------------
    Q(:,:,i) = (R + G);
    
end

% basis for expansion of Q operator
%--------------------------------------------------------------------------
b     = spm_dctmtx(nx,fix(nx/2));
[Db,Dx]  = spm_dctmtx(nx,fix(nx/2),'diff',dx);
for i = 1:n
    bdb{i,1} = b'*Dx{i};
end
Db    = full(Db);


for i = 1:n
    for j = 1:n
        Qb(i,j,:) = b'*reshape(Q(i,j,:),nX,1);
    end
end


dQd0  = cell(n,n); [dQd0{:}] = deal(sparse(nX,nX));
dQdB  = {};
for i = 1:n
    for j = (i + 1):n
        for k = 1:size(b,2)
            
            dQdb      = dQd0;
            dQdb{i,j} =  spd(b(:,k));
            dQdb{j,i} = -spd(b(:,k));
            dQdB{end + 1}    = spm_cat(dQdb);
            
        end
    end
end

nb    = size(b,2);
nB    = numel(dQdB);
dqdb  = {};
for i = 1:n
    for j = (i + 1):n
        for k = 1:size(b,2)
            
            dqdb{end + 1}    = spm_zeros(Qb);
            dqdb{end}(i,j,k) =  1;
            dqdb{end}(j,i,k) = -1;
            
        end
    end
end

dLdb  = cell(1,nB); [dLdb{:}] = deal(sparse(nb,n));
for i = 1:n
    for j = 1:n
        for k = 1:nB
            dLdb{k}(:,i) = dLdb{k}(:,i) - bdb{j}*reshape(dqdb{k}(i,j,:),nb,1);
        end
    end
end

   
dEdb  = zeros(nX*n,nB);
QQ    = cell(n,n);
for q = 1:8
     

    % correction terms state-dependent solenoidal flow
    %----------------------------------------------------------------------
    L     = zeros(nb,n);
    for i = 1:n
        for j = 1:n
            L(:,i) = L(:,i) - bdb{j}*reshape(Qb(i,j,:),nb,1);
        end
    end
    
    % Kroneckor form of Q
    %----------------------------------------------------------------------
    for i = 1:n
        for j = 1:n
            QQ{i,j} = spd(b*squeeze(Qb(i,j,:)));
        end
    end
    QZ  = spm_cat(QQ);
        
    % potential gradients f = QZ*DZ*b*b'*q0 - L
    %----------------------------------------------------------------------
    Z   = QZ*Db;
    Y   = spm_vec(f' + b*L);
    q0  = (Z'*Z + speye(nb,nb)*exp(-16))\(Z'*Y);

    
    % dEdb  = dLdb - dQdb*DZ*q0;
    %----------------------------------------------------------------------
    Dbq0 = Db*q0;
    for i = 1:nB
        dEdb(:,i) = dQdB{i}*Dbq0 - spm_vec(b*dLdb{i});
    end
    
    % dEdb  = dLdb - dQdb*DZ*q0;
    %----------------------------------------------------------------------
    E     = Y - Z*q0;
    B     = (dEdb'*dEdb + speye(nB,nB)*exp(-16))\(dEdb'*E);
    for i = 1:numel(dqdb)
        Qb = Qb + dqdb{i}*B(i);
    end
    
    subplot(3,2,1)
    plot(spm_vec(f'),Z*q0 - spm_vec(b*L),'.','MarkerSize',1),drawnow
    disp(norm(E(:)))
    subplot(3,2,2)
    p0  = reshape(b*q0,nx); imagesc(squeeze(sum(p0,3)))
    
end


for i = 1:n
    for j = 1:n
        Q(i,j,:) = b*reshape(Qb(i,j,:),nb,1);
    end
end

    
% precision (inverse covariance) of steady-state density
%--------------------------------------------------------------------------
H     = zeros(n,n,nX);
E     = zeros(n,nX);
for i = 1:nX
    
    % precision (inverse covariance) of steady-state density
    %----------------------------------------------------------------------
    H(:,:,i) = -full(Q(:,:,i)\J(:,:,i));
    
    % eigenvalues of Jacobian
    %----------------------------------------------------------------------
    E(:,i)   = sort(eig(J(:,:,i)),'descend','ComparisonMethod','real');
    
end

% expected Euclidian norm of Hessian (first state is a blanket state)
%--------------------------------------------------------------------------
p0  = spm_softmax(b*q0);
H   = spm_dot(H.^2,p0);
J   = spm_dot(J.^2,p0);
E   = spm_dot(E   ,p0);

% 2+abs(E(1)+E(2))/abs(E(3))

% reshape nonequilibrium steady-state density
%--------------------------------------------------------------------------
p0  = reshape(p0,nx);

% predicted flow     Y   = spm_vec((f*b)' + L);
%--------------------------------------------------------------------------
F   = reshape(Z*q0,nX,n) - b*L;


return


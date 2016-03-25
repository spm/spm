function [K0,K1,K2,H1] = spm_kernels(varargin)
% returns global Volterra kernels for a MIMO Bilinear system
% FORMAT [K0,K1,K2] = spm_kernels(M,P,N,dt)            - output kernels
% FORMAT [K0,K1,K2] = spm_kernels(M0,M1,N,dt)          - state  kernels
% FORMAT [K0,K1,K2] = spm_kernels(M0,M1,L1,N,dt)       - output kernels (1st)
% FORMAT [K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt)    - output kernels (2nd)
%
% M,P   - model structure and parameters; 
%         or its bilinear reduction:
%
% M0    - (n x n)     df(q(0),0)/dq                    - n states
% M1    - {m}(n x n)  d2f(q(0),0)/dqdu                 - m inputs
% L1    - (l x n)     dldq                             - l outputs
% L2    - {m}(n x n)  dl2dqq
%
% N     - kernel depth       {intervals}
% dt    - interval           {seconds}
%
% Volterra kernels:
%---------------------------------------------------------------------------
% K0    - (1 x l)             = K0(t)         = y(t)
% K1    - (N x l x m)         = K1i(t,s1)     = dy(t)/dui(t - s1)
% K2    - (N x N x l x m x m) = K2ij(t,s1,s2) = d2y(t)/dui(t - s1)duj(t - s2)
%
%___________________________________________________________________________
% Returns Volterra kernels for bilinear systems of the form
%
%         dq/dt   = f(q,u) = M0*q + M1{1}*q*u1 + ... M1{m}*q*um
%            y(i) = L1(i,:)*q + q'*L2{i}*q
%
% where q = [1 x(t)] are the states augmented with a constant term
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_kernels.m 6755 2016-03-25 09:48:34Z karl $
 
 
% assign inputs
%--------------------------------------------------------------------------
if nargin == 4
 
    M0 = varargin{1};
    M1 = varargin{2};
    N  = varargin{3};
    dt = varargin{4};
 
elseif nargin == 5
 
    M0 = varargin{1};
    M1 = varargin{2};
    L1 = varargin{3};
    N  = varargin{4};
    dt = varargin{5};
 
elseif nargin == 6
 
    M0 = varargin{1};
    M1 = varargin{2};
    L1 = varargin{3};
    L2 = varargin{4};
    N  = varargin{5};
    dt = varargin{6};
end
 
% bilinear reduction if necessary
%--------------------------------------------------------------------------
if isstruct(M0)
    [M0,M1,L1,L2] = spm_bireduce(M0,M1);
end
 
 
% Volterra kernels for bilinear systems
%==========================================================================
 
% make states the outputs (i.e. remove constant) if L1 is not specified
%--------------------------------------------------------------------------
try
    L1;
catch
    L1 = speye(size(M0));
    L1 = L1(2:end,:);
end
try
    L2;
catch
    L2 = [];
end
 
% parameters
%--------------------------------------------------------------------------
N     = fix(N);                     % kernel depth
n     = size(M0,1);                 % state variables
m     = size(M1,1);                 % inputs
l     = size(L1,1);                 % outputs
H1    = zeros(N,n,m);
K1    = zeros(N,l,m);
K2    = zeros(N,N,l,m,m);
M0    = full(M0);
 
% pre-compute matrix exponentials
%--------------------------------------------------------------------------
e1    = expm( dt*M0);
e2    = expm(-dt*M0);
for p = 1:m
    M{1,p} = e1*M1{p}*e2;
end
ei    = e1;
for i = 2:N
    ei    = e1*ei;
    for p = 1:m
        M{i,p} = e1*M{i - 1,p}*e2;
    end
end
 
 
% check for convergence
%--------------------------------------------------------------------------
q     = 0;
for p = 1:m
    q = q | norm(M{N,p},'inf') > norm(M{1,p},'inf');
    q = q | norm(M{N,p},'inf') > exp(16);
    q = q | isnan(norm(M{N,p}));
end
 
% if necessary, remove unstable modes and repeat with a pseudoinverse
%--------------------------------------------------------------------------
if q
    M0    = spm_bilinear_condition(M0);
    e1    = expm(dt*M0);
    ei    = 1;
    for i = 1:N
        ei    = e1*ei;
        ie    = spm_pinv(ei);
        for p = 1:m
            M{i,p} = ei*M1{p}*ie;
        end
    end
end
 
 
% 0th order kernel
%--------------------------------------------------------------------------
X0    = sparse(1,1,1,n,1);
if nargout > 0
    H0    = ei*X0;
    K0    = L1*H0;
end
 
 
% 1st order kernel
%--------------------------------------------------------------------------
if nargout > 1
    for p = 1:m
        for i = 1:N
            H1(i,:,p) = M{i,p}*H0;
            K1(i,:,p) = H1(i,:,p)*L1';
        end
    end
end
 
 
% 2nd order kernels
%--------------------------------------------------------------------------
if nargout > 2
    for p = 1:m
        for q = 1:m
            for j = 1:N
                H  = L1*M{j,q}*H1((j:N),:,p)';
                K2(j,j:N,:,q,p) = H';
                K2(j:N,j,:,p,q) = H';
 
            end
        end
    end
 
    if isempty(L2), return, end
 
    % add output nonlinearity
    %----------------------------------------------------------------------
    for i = 1:m
        for j = 1:m
            for p = 1:l
                K2(:,:,p,i,j) = K2(:,:,p,i,j) + H1(:,:,i)*L2{p}*H1(:,:,j)';
            end
        end
    end
end

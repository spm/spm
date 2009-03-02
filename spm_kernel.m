function [K1,K0,K2] = spm_kernel(M0,P,N,dt)
% returns the first Volterra kernel for a MIMO Bilinear system
% FORMAT [K1,K0,K2] = spm_kernel(M,P,N,dt)
%
% M     - model structure
% P     - model parameters
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
% see also spm_kernels
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_kernel.m 2804 2009-03-02 12:03:00Z karl $


% bilinear reduction if necessary
%--------------------------------------------------------------------------
if nargout > 2
    [M0,M1,L1,L2] = spm_bireduce(M0,P);
else
    [M0,M1,L1]    = spm_bireduce(M0,P);
end


% Volterra kernels for bilinear systems
%==========================================================================

% make states the outputs (i.e. remove constant) if L1 is not specified
%--------------------------------------------------------------------------
try
    L2;
catch
    L2 = [];
end

% parameters
%--------------------------------------------------------------------------
N     = fix(N);                     % kernel depth
n     = size(M0,1);                 % state variables
m     = size(M1,2);                 % inputs
l     = size(L1,1);                 % outputs
H1    = zeros(N,n,m);
K1    = zeros(N,l,m);
K2    = zeros(N,N,l,m,m);
M0    = full(M0);


% pre-compute exponentials
%--------------------------------------------------------------------------
e1    = sparse(expm( dt*M0));
e2    = sparse(expm(-dt*M0));
for p = 1:m
    M{1,p} = e1*M1{p}*e2;
end
for i = 2:N
    for p = 1:m
        M{i,p} = e1*M{i - 1,p}*e2;
    end
end

% 0th order kernel
%--------------------------------------------------------------------------
X0    = sparse(1,1,1,n,1);
if nargout > 0
    H0    = e1^N*X0;
    K0    = L1*H0;
end


% 1st order kernel
%--------------------------------------------------------------------------
for p = 1:m
    for i = 1:N

        % 1st order kernel
        %------------------------------------------------------------------
        H1(i,:,p) = M{i,p}*H0;
        K1(i,:,p) = H1(i,:,p)*L1';
    end
end


% 2nd order kernels
%--------------------------------------------------------------------------
if nargout > 2
    for p = 1:m
        for q = 1:m
            for j = 1:N

                % 2nd order kernel
                %----------------------------------------------------------
                H  = L1*M{j,q}*H1([j:N],:,p)';
                K2(j,[j:N],:,q,p) = H';
                K2([j:N],j,:,p,q) = H';

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

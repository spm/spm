function [J,z,v,s] = spm_A_reduce(J,x,T,N)
% Reduction of Markovian partition
% FORMAT [J,z,v,s] = spm_A_reduce(J,x,T,N)
% J  - Jacobian (x)
% x  - {3 x n}  particular partition of states
% T  - eigenvalue threshold to retain eigenvectors [default: 8]
% N  - maximum number to retain [default: 8]
%
% J  - Jacobian (z)
% z  - {1 x n} partition of states at the next level
% v  - {1 x n} eigenvector (adiabatic) operator
% s  - {1 x n} eigenvalues
%
% Adiabatic reduction operator (R)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
nx    = size(x,2);                  % number of partitions
if nargin < 3
    T = 8;                          % adiabatic threshold (Hz)
end
if nargin < 4
    N = 8;                          % maximum number
end


% reduction
%--------------------------------------------------------------------------
for i = 1:nx
    
    % Lyapunov exponents (eigensolution) for this partition
    %----------------------------------------------------------------------
    y{i}  = spm_vec(x(1:2,i));
    Jii   = full(J(y{i},y{i}));
    [e,r] = eig(Jii);
    r     = diag(r);
    [d,j] = sort(real(r),'descend');
    
    % Adiabatic threshold
    %----------------------------------------------------------------------
    t     = max(-T,d(min(end,N + 1)));
    n(i)  = sum(d > t);
    s{i}  = r(  j(1:n(i)));
    v{i}  = e(:,j(1:n(i)));
    u{i}  = pinv(v{i});
    
end
for i = 1:nx
    for j = 1:nx
        Jij    = full(J(spm_vec(x(1:2,i)),spm_vec(x(1:2,j))));
        A{i,j} = u{i}*Jij*v{j};
    end
    z{i}       = sum(n(1:(i - 1))) + (1:n(i));
end
J  = spm_cat(A);

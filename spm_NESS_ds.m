function [dS,G,Q,L] = spm_NESS_ds(Sp,P,x)
% Generate changes in log density (coefficients or at x)
% FORMAT [dS,G,Q,L] = spm_NESS_ds(Sp,P)
% FORMAT [ds,G,Q,L] = spm_NESS_ds(Sp,P,x)
%--------------------------------------------------------------------------
% Sp      - polynomial coefficients of initial potential
% P.Qp    - polynomial coefficients of solenoidal operator
% P.Sp    - polynomial coefficients of final potential
% P.G     - amplitude of random fluctuations
%
% x       - sample points and state space
%
% dS      - time derivative of polynomial coefficients of potential
% ds      - time derivative of potential at x
% G       - dissipation operator
% Q       - solenoidal operator
% L       - correction term for derivatives of solenoidal flow
%
% This routine assumes that K = 3; i.e., the log density is second order in
% the states (Laplace assumption). if call with two arguments the time
% derivatives of the (second-order) polynomial coefficients of the log
% density are returned. If called with three arguments, the time derivative
% of the log density at the specified points in state space are returned.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


% get sample points
%==========================================================================
if nargin < 3
    n    = size(P.G,1);
    N    = 3;
    for i = 1:n
        x{i} = linspace(-1,1,N);
    end
end

% assume the log density is second-order in the states (Laplace assumption)
%--------------------------------------------------------------------------
[b,D,H] = spm_polymtx(x,3);
[nX,nb] = size(b);
n       = numel(x);

% derivatives of flow operator Q
%--------------------------------------------------------------------------
Q     = zeros(nX,n,n,'like',P.Qp);
Qp    = zeros(nb,n,n,'like',P.Qp);
k     = 0;
for i = 1:n
    for j = i:n
        k         = k + 1;
        Qp(:,i,j) = P.Qp((1:nb) + (k - 1)*nb);
        
        % diagonal term
        %------------------------------------------------------------------
        Q(:,i,j) = b*Qp(:,i,j);
        
        % skew symmetric terms
        %--------------------------------------------------------------
        if i < j
            Q(:,j,i) = -Q(:,i,j);
        end
    end
end


% correction term L
%--------------------------------------------------------------------------
L     = zeros(nX,n,'like',P.Qp);
for i = 1:n
    for j = 1:n
        L(:,i) = L(:,i) - D{j}*Qp(:,i,j);
    end
end

% dS potential difference
%--------------------------------------------------------------------------
G     = P.G;
ds    = zeros(nX,1,'like',b);
for i = 1:n
    
    % correction term
    %----------------------------------------------------------------------
    ds = ds - L(:,i).*D{i}*(P.Sp - Sp);
    
    % curvature term
    %----------------------------------------------------------------------
    ds = ds + G(i,i)*H{i,i}*(P.Sp - Sp);
    
    % curvature term
    %----------------------------------------------------------------------
    ds = ds - G(i,i)*(D{i}*Sp).*(D{i}*(P.Sp - Sp));
    
    % solenoidal term
    %----------------------------------------------------------------------
    for j = 1:n
        ds = ds + (D{i}*P.Sp).*Q(:,i,j).*(D{j}*Sp);
    end
end

% return time derivatives of coefficients or log density
%--------------------------------------------------------------------------
if nargin < 3
    dS = b\ds;
else
    dS = ds;
end

return



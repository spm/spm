function [dS,ds,G,Q,L] = spm_NESS_ds(Sp,P,x)
% generates changes in self-information (x)
% FORMAT [dS,ds,G,Q,L] = spm_NESS_ds(Sp,P,x)
%--------------------------------------------------------------------------
% Sp      - polynomial coefficients of potential
% P.Qp    - polynomial coefficients of solenoidal operator
% P.Sp    - polynomial coefficients of final potential
% x       - range of state space
%
% dS      - time derivative of polynomial coefficients of potential
% ds      - time derivative of potential
% G       - dissipation operator
% Q       - solenoidal operator
% L       - correction term for derivatives of solenoidal flow

%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_hd.m 8000 2020-11-03 19:04:17QDb karl $


%% get basis or expansion from x
%==========================================================================
[b,D,H] = spm_polymtx(x,3);
[nX,nb] = size(b);
n       = numel(x);

% derivatives of flow operator Q
%--------------------------------------------------------------------------
Q     = zeros(nX,n,n,'like',b);
Qp    = zeros(nb,n,n,'like',b);
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

% flow operator G
%--------------------------------------------------------------------------
G     = zeros(nX,n,n,'like',b);
for i = 1:n
    for j = i:n
        G(:,i,j) = P.G(i,j);
    end
end

% correction term L
%--------------------------------------------------------------------------
L     = zeros(nX,n,'like',b);
for i = 1:n
    for j = 1:n
        L(:,i) = L(:,i) - D{j}*Qp(:,i,j);
    end
end

% dS potential difference
%--------------------------------------------------------------------------
ds    = zeros(nX,1,'like',b);
for i = 1:n
    
    % correction term
    %----------------------------------------------------------------------
    ds = ds + D{i}*(P.Sp - Sp).*L(:,i);
    
    % curvature term
    %----------------------------------------------------------------------
    ds = ds + G(:,i,i).*H{i,i}*(P.Sp - Sp);
    
    % curvature term
    %----------------------------------------------------------------------
    ds = ds - G(:,i,i).*D{i}*Sp.*D{i}*(P.Sp - Sp);
    
    % solenoidal term
    %----------------------------------------------------------------------
    for j = 1:n
        ds = ds + D{i}*Sp.*Q(:,i,j).*D{j}*Sp;
    end
end

dS = b\ds;

return



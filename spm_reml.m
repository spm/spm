function [Ce,h,W,u,F] = spm_reml(Cy,X,Q);
% ReML estimation of covariance components from Cov{y}
% FORMAT [Ce,h,W,u,F] = spm_reml(Cy,X,Q);
%
% Cy  - (m x m) data covariance matrix y*y'/n  {y = (m x n) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
%
% Ce  - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) hyperparameters
% W   - (q x q) W*n = precision of hyperparameter estimates
% u   - {1 x p} estimable components C{i} = u(1,i)*Q{1} + u(2,i)*Q{2} +...
%
% F   - [-ve] free energy F = log evidence = p(y|X,Q) = ReML objective
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% esitmates.  An SVD of the curvature is used to realign the covariance
% components with its principal axes to ensure all components are estimated
% with finite precision.
%_______________________________________________________________________
% @(#)spm_reml.m	2.33 John Ashburner, Karl Friston 04/09/27

% Tolerances 
%-----------------------------------------------------------------------
TOL   = 1e-6;       % for convergence [norm of gradient df/dh]
TOS   = 1e-6;	    % for SVD of curvature [ddf/dhdh]

% ensure X is not rank deficient
%-----------------------------------------------------------------------
X     = full(X);
X     = orth(X);
X     = sparse(X);

% find estimable components (encoded in the precision matrix W)
%-----------------------------------------------------------------------
m     = length(Q);
n     = length(Cy);
W     = zeros(m,m);
for i = 1:m
    RQ{i}  = Q{i} - X*(X'*Q{i});
end
for i = 1:m
    for j = i:m
        dFdhh  = sum(sum(RQ{i}.*RQ{j}'));
        W(i,j) = dFdhh;
        W(j,i) = dFdhh;
    end
end

% eliminate inestimable components
%-----------------------------------------------------------------------
[u s] = spm_svd(W,TOS);
for i = 1:size(u,2)
    C{i}  = sparse(n,n);
    for j = 1:m
        C{i} = C{i} + Q{j}*u(j,i);
    end
end
Q     = C;

% initialize hyperparameters
%-----------------------------------------------------------------------
m     = length(Q);
dFdh  = zeros(m,1);
W     = zeros(m,m); 
C     = [];
for i = 1:m
    C = [C Q{i}(:)];
end
R     = speye(n,n) - X*X';
R     = R*Cy*R';
h     = inv(C'*C)*(C'*R(:));
dh    = sparse(m,1);

% and Ce		
%-----------------------------------------------------------------------
Ce    = sparse(n,n);
for i = 1:m
    Ce = Ce + h(i)*Q{i};
end

% Iterative EM
%-----------------------------------------------------------------------
for k = 1:32
    
    % Q are variance components		
    %-------------------------------------------------------------------
    dC    = sparse(n,n);
    for i = 1:m
        dC = dC + dh(i)*Q{i};
    end
    
    % Check Ce is positive semi-definite
    %-------------------------------------------------------------------
    if any(diag(Ce + dC) < 0)
        Ce = Ce + dC/2;
        h  = h  + dh/2;
    else
        Ce = Ce + dC;
        h  = h  + dh;
    end
    iCe   = inv(Ce);
    
    
    % E-step: conditional covariance cov(B|y) {Cby}
    %===================================================================
    iCeX  = iCe*X;
    Cby   = inv(X'*iCeX);
    
    % M-step: ReML estimate of hyperparameters 
    %===================================================================
    
    % Gradient dFd/h (first derivatives)
    %-------------------------------------------------------------------
    P     = iCe  - iCeX*Cby*iCeX';
    PCy   = Cy*P'- speye(n,n);
    for i = 1:m
        
        % dF/dh = -trace(dF/diCe*iCe*Q{i}*iCe) = 
        %---------------------------------------------------------------
        PQ{i}   = P*Q{i};
        dFdh(i) = sum(sum(PCy.*PQ{i}))/2;
    end
    
    % Expected curvature E{ddF/dhh} (second derivatives)
    %-------------------------------------------------------------------
    for i = 1:m
        for j = i:m
            
            % ddF/dhh = -trace{P*Q{i}*P*Q{j}}
            %-----------------------------------------------------------
            dFdhh  = sum(sum(PQ{i}.*PQ{j}))/2;
            W(i,j) = dFdhh;
            W(j,i) = dFdhh;
        end
    end
    
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %-------------------------------------------------------------------
    dh    = pinv(W)*dFdh(:);
    
    % Convergence
    %===================================================================
    w     = (dh'*dh)/(h'*h);
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(w));
    if w < TOL, break, end

end

% project hyperparameter esimates and precision back to user-specifed Q
%-----------------------------------------------------------------------
if nargout > 1
    h     = u*h;
    W     = u*W*u';
end

% log evidence = ln p(y|X,Q) = ReML objective = F
%-----------------------------------------------------------------------
if nargout > 3
     F = -trace(P'*Ce'*P*Cy)/2 ...
         -n*log(2*pi)/2 ...
         -spm_logdet(Ce)/2 ...
         +spm_logdet(Cby)/2;
end

function [C,h,Ph,F] = spm_reml(YY,X,Q,N,OPT);
% ReML estimation of covariance components from y*y'
% FORMAT [C,h,Ph,F] = spm_reml(YY,X,Q,N,[OPT]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x n) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% OPT = 1 : log-normal hyper-parameterisation
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of h [or log(h), if OPT(1)]
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% John Ashburner & Karl Friston
% $Id: spm_reml.m 249 2005-10-05 17:58:37Z karl $
 
% assume a single sample if not specified
%--------------------------------------------------------------------------
try
    N;
catch
    N  = 1;
end
 
% assume OPT = [0 0]
%--------------------------------------------------------------------------
try
    OPT;
catch
    OPT = 0;
end
 
% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(length(Q{1}),1);
end
X     = orth(full(X));
[n p] = size(X);
 
% initialise h
%--------------------------------------------------------------------------
m     = length(Q);
dh    = sparse(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);

% log-transform if necessary
%--------------------------------------------------------------------------
if OPT
    [C,h] = spm_reml(YY,X,Q,N);
    fprintf('%s:\n','Applying log-normal hyperpriors');
    h     = log(max(h,1e-6));

else
    h     = ones(m,1);
end
 
 
% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:32
    
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        if OPT
            C = C + Q{i}*exp(h(i));
        else
            C = C + Q{i}*h(i);
        end
    end
    iC    = inv(C);
 
    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iCX   = iC*X;
    Cq    = pinv(X'*iCX);
    XCXiC = X*Cq*iCX';
 
    % M-step: ReML estimate of hyperparameters
    %======================================================================
 
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    P     = iC - iC*XCXiC;
    U     = speye(n) - P*YY/N;
    for i = 1:m
 
        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        PQ{i}     = P*Q{i};
        if OPT
            PQ{i} = PQ{i}*exp(h(i));
        end
        dFdh(i)   = -sum(sum(PQ{i}.*U'))*N/2;

    end
 
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m
 
            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j) = -sum(sum(PQ{i}.*PQ{j}'))*N/2;
            dFdhh(j,i) =  dFdhh(i,j);
 
        end
    end
    
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    Ph    = -dFdhh;
    dh    = -pinv(dFdhh)*dFdh;
    h     = h + dh;
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    w     = dFdh'*dh;
    if w < 1e-2, break, end
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(w));

 
end
 
% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
if nargout > 3
 
    F = - trace(C*P*YY*P)/2 ...
        - N*n*log(2*pi)/2 ...
        - N*spm_logdet(C)/2 ...
        + N*spm_logdet(Cq)/2 ...
        -   spm_logdet(Ph)/2 ...
        - N*p/2 - m/2;
end
 
% return exp(h) if log-normal hyperpriors
%--------------------------------------------------------------------------
if OPT
    h  = exp(h);
end

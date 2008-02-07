function [C,h,Ph,F] = spm_fn_reml(YY,X,Q,N,hE,K);
% ReML estimation of covariance components from y*y'
% FORMAT [C,h,Ph,F] = spm_fn_reml(YY,X,Q,N,hE,K);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - inline function or script C = Q(h,m)
% N   - number of samples
%
% hE  - prior expectation (& starting esitmate for Q(h,m))
% K   - maxmium number of iterations
%
% C   - (m x m) estimated errors: C = Q(h)
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of h
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Karl Friston
% $Id: spm_fn_reml.m 1143 2008-02-07 19:33:33Z spm $

% assume a single sample if not specified
%--------------------------------------------------------------------------
try
    N;
catch
    N  = 1;
end

% default number of iterations
%--------------------------------------------------------------------------
try
    K;
catch
    K  = 64;
end


% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(length(YY),1);
else
    X = orth(full(X));
end

% initialise h
%--------------------------------------------------------------------------
n     = length(YY);
m     = length(hE);
h     = hE;
dh    = zeros(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
L     = zeros(m,m);


% specify hyperpriors
%--------------------------------------------------------------------------
hP    = speye(m,m)/exp(32);



% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:K

    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = feval(Q,h,n);
    iC    = inv(C + speye(n,n)/exp(32));

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
        PQ{i}     = P*spm_diff(Q,h,n,1);
        dFdh(i)   = -trace(PQ{i}*U)*N/2;

    end

    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m

            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j) = -trace(PQ{i}*PQ{j})*N/2;
            dFdhh(j,i) =  dFdhh(i,j);

        end
    end
    
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;
    
    % update regulariser
    %----------------------------------------------------------------------
    L  = speye(m,m)*norm(dFdhh,1)/128;
    
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    Ph    = -dFdhh;
    dh    = -inv(dFdhh - L)*dFdh;

    % preclude numerical overflow
    %----------------------------------------------------------------------
    h     = h + dh;
    if nargin > 4
        h = min(h, 32);
        h = max(h,-32);
    end
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    dF    = dFdh'*dh;
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(dF));
    if dF < 1e-1, break, end

end

% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
if nargout > 3
    
    F = - trace(C*P*YY*P)/2 ...
        - e'*hP*e/2 ...
        - N*n*log(2*pi)/2 ...
        - N*spm_logdet(C)/2 ...
        + N*spm_logdet(Cq)/2 ...
        -   spm_logdet(Ph)/2 ...
        +   spm_logdet(hP)/2;
end

    

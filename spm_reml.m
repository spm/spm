function [C,h,Ph,F] = spm_reml(YY,X,Q,N,hE,hC);
% ReML estimation of covariance components from y*y'
% FORMAT [C,h,Ph,F] = spm_reml(YY,X,Q,N,[hE,hC]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% hE  - hyperprior expectation in log-space [default = 0]
% hC  - hyperprior covariance  in log-space [default = 4]
%     - If specified hE induces log-normal hyperpriors
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of h [or log(h), if hE(1)]
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner & Karl Friston
% $Id: spm_reml.m 666 2006-10-25 14:05:00Z karl $

% assume a single sample if not specified
%--------------------------------------------------------------------------
try, N; catch, N  = 1;  end

% default number of iterations
%--------------------------------------------------------------------------
try, K; catch, K  = 64; end

% hyperpriors; if not specified
%--------------------------------------------------------------------------
try, hE;            catch, hE = 0;   end
try, hP = inv(hC);  catch, hP = 1/4; end

% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(length(Q{1}),1);
else
    X = orth(full(X));
end

% initialise h
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
h     = zeros(m,1);
dh    = zeros(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
L     = zeros(m,m);


% initialise and specify hyperpriors
%--------------------------------------------------------------------------
if nargin > 4
    LY    = log(normest(YY)) - log(N);
    for i = 1:m
        h(i) = LY - log(normest(Q{i}));
    end
else
    for i = 1:m
        h(i) = any(diag(Q{i}));
    end
    hP  = speye(m,m)/exp(32);
end


% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:K

    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        if nargin > 4
            C = C + Q{i}*exp(h(i));
        else
            C = C + Q{i}*h(i);
        end
    end
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
        PQ{i}     = P*Q{i};
        if nargin > 4
            PQ{i} = PQ{i}*exp(h(i));
        end
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
    if ~rem(k,8)
       L  = speye(m,m)*norm(dFdhh)/128;
    end
    
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    Ph    = -dFdhh;
    dh    = -pinv(dFdhh - L)*dFdh;

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

% return exp(h) if log-normal hyperpriors
%--------------------------------------------------------------------------
if nargin > 4
    h  = exp(h);
end
    

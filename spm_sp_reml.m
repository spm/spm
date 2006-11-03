function [C,h,Ph,F] = spm_sp_reml(YY,X,Q,N,hE,hC);
% ReML estimation of covariance components from y*y'
% FORMAT [C,h,Ph,F] = spm_sp_reml(YY,X,Q,N,[hE,hC]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components Q.q = eigenvectors; Q.v = eigen
%               values
% N   - number of samples
%
% hE  - hyperprior expectation in log-space [default = 0]
% hC  - hyperprior covariance  in log-space [default = 4]
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
% $Id: spm_reml.m 615 2006-09-08 16:16:06Z karl $

% assume a single sample if not specified
%--------------------------------------------------------------------------
try, N; catch, N  = 1;  end

% default number of iterations
%--------------------------------------------------------------------------
try, K; catch, K  = 256; end

% hyperpriors; if not specified
%--------------------------------------------------------------------------
try, hE;            catch, hE = 0;   end
try, hP = inv(hC);  catch, hP = 1/4; end

% call spm_reml as default
%--------------------------------------------------------------------------
try
    [C,h,Ph,F] = spm_reml(YY,X,Q,N,hE,hC);
    return
end

% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(length(Q{1}),1);
else
    X = orth(full(X));
end

% find bases of Q if necessary
%--------------------------------------------------------------------------
for i = 1:length(Q)
    if isstruct(Q{i})
        try
            if ~isfield(Q{i},'v');
                Q{i}.v = ones(size(Q{i}.q,2),1);
            end
        catch
            warndlg('Please specify components with Q.q (basis) and Q.v');
            return
        end
    else
        [q v] = spm_svd(Q{i},exp(-16));
        C.q   = q;
        C.v   = diag(v);
        Q{i}  = C;
    end
end

% compute basis and dsdh
%--------------------------------------------------------------------------
m     = length(Q);
q     = cell(1,m);
v     = cell(1,m);
for i = 1:m
    q{i} = Q{i}.q;
    v{i} = Q{i}.v(:);
end
q     = spm_cat(q);
dedh  = spm_cat(diag(v));

% initialise h
%--------------------------------------------------------------------------
[n s] = size(q);
h     = zeros(m,1);
dh    = zeros(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
L     = zeros(m,m);

% initialise hyperparameters and specify hyperpriors
%--------------------------------------------------------------------------
h     = sparse(m,1);

% pre-comooute bases
%--------------------------------------------------------------------------
qq    = cell(s,1);
for i = 1:s
    qq{i} = q(:,i)*q(:,i)';
end

% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:K

    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    e     = dedh*exp(h);
    for i = 1:s
        C = C + qq{i}*e(i);
    end
    iC    = inv(C + speye(n,n)/exp(32));

    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iCX   = iC*X;
    Cq    = inv(X'*iCX);
    XCXiC = X*Cq*iCX';

    % M-step: ReML estimate of hyperparameters
    %======================================================================

    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    U     = iC - iC*XCXiC;
    W     = U*(YY/N - C)*U';
    P     = q'*U*q;

    % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
    %----------------------------------------------------------------------
    dFde  =  N/2*sum(q.*(W*q))';
    dFdee = -N/2*P.*P';

    % dF/dhh = -trace{P*Q{i}*P*Q{j}}
    %----------------------------------------------------------------------
    dhdh  = dedh*diag(exp(h));
    dFdh  = dhdh'*dFde;
    dFdhh = dhdh'*dFdee*dhdh;
    
    % add hyperpriors
    %----------------------------------------------------------------------
    dFdh  = dFdh  - hP*(h - hE);
    dFdhh = dFdhh - hP;
    
    % update regulariser
    %----------------------------------------------------------------------
    L     = speye(m,m)*norm(dFdhh,1)/256;
 
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    Ph    = -dFdhh;
    dh    = -inv(dFdhh - L)*dFdh;

    % preclude numerical overflow
    %----------------------------------------------------------------------
    h     = h + dh;
    h     = min(h, 32);
    h     = max(h,-32);
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    dF    = dFdh'*dh;
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(dF));
    if dF < 1e-1, break, end

end

% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
if nargout > 3
    
    F = - trace(C*U*YY*U')/2 ...
        - (h - hE)'*hP*(h - hE)/2 ...
        - N*n*log(2*pi)/2 ...
        - N*spm_logdet(C)/2 ...
        + N*spm_logdet(Cq)/2 ...
        -   spm_logdet(Ph)/2 ...
        +   spm_logdet(hP)/2;
end

% return exp(h) if log-normal hyperpriors
%--------------------------------------------------------------------------
h  = exp(h);

    

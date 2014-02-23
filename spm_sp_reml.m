function [C,h,Ph,F,Fa,Fc] = spm_sp_reml(YY,X,Q,N,hE);
% ReML estimation of covariance components from y*y' (for sparse patterns)
% FORMAT [C,h,Ph,F,Fa,Fc] = spm_sp_reml(YY,X,Q,N);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} components Q.q = eigenvectors; Q.v = eigenvalues
%               or (m x n) matrix of n basis functions
% N   - number of samples
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of log(h)
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates,. using uninformative hyperpriors (this is effectively an ARD
% scheme).  The specification of components differs from spm_reml and
% spm_reml_sc.
%
%__________________________________________________________________________
%
% SPM ReML routines:
%
%      spm_reml:    no positivity constraints on covariance parameters
%      spm_reml_sc: positivity constraints on covariance parameters
%      spm_sp_reml: for sparse patterns (c.f., ARD)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_sp_reml.m 5892 2014-02-23 11:00:16Z karl $
 
% assume a single sample if not specified
%--------------------------------------------------------------------------
try, N; catch, N = 1;   end
 
% default number of iterations
%--------------------------------------------------------------------------
try, K; catch, K = 128; end
 
% number of hyperparameters
%--------------------------------------------------------------------------
n     = length(YY);
m     = length(Q);
h     = zeros(m,1);
dh    = zeros(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
 
% uninformative hyperpriors
%--------------------------------------------------------------------------
hE    = sparse(m,1) - 32;
hC    = speye(m,m)*256;
hP    = inv(hC);

% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
    R = speye(n,n);
else
    X = orth(full(X));
    R = speye(n,n) - X*X';
end

% convert bssis to struct if necessary
%--------------------------------------------------------------------------
if ~iscell(Q)
    q     = cell(1,m);
    for i = 1:m
        q{i}.q = spm_en(Q(:,i));
    end
    Q     = q;
end
    
% find bases of Q if necessary
%--------------------------------------------------------------------------
for i = 1:m
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
        [q,v] = spm_svd(Q{i});
        C.q   = q;
        C.v   = diag(v);
        Q{i}  = C;
    end
end
 
% scale YY
%--------------------------------------------------------------------------
sY    = n*trace(YY)/N;
YY    = YY/sY;
 
% scale Q
%--------------------------------------------------------------------------
for i = 1:m
    sh(i,1) = n*trace(Q{i}.q'*R*Q{i}.q*diag(Q{i}.v));
    Q{i}.q  = Q{i}.q/sqrt(sh(i));
end
 
% compute basis and dsdh
%--------------------------------------------------------------------------
q     = cell(1,m);
v     = cell(1,m);
for i = 1:m
    q{i} = Q{i}.q;
    v{i} = Q{i}.v(:);
end
q     = spm_cat(q);
dedh  = spm_cat(spm_diag(v));
 
 
% pre-compute bases
%--------------------------------------------------------------------------
[n s] = size(q);
qq    = cell(s,1);
for i = 1:s
    qq{i} = q(:,i)*q(:,i)';
end
 
% ReML (EM/VB)
%--------------------------------------------------------------------------
dF    = Inf;
as    = 1:m;
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
    iCX    = iC*X;
    if ~isempty(X)
        Cq = inv(X'*iCX);
    else
        Cq = sparse(0);
    end
 
    % M-step: ReML estimate of hyperparameters
    %======================================================================
 
    % select relevant patterns
    %----------------------------------------------------------------------
    rel   = any(dedh(:,as),2);
    
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    U     = iC - iCX*Cq*iCX';
    W     = U*(YY/N - C)*U';
    qr    = q(:,rel);
    P     = qr'*U*qr;
 
    % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
    %----------------------------------------------------------------------
    dFde  =  N/2*sum(qr.*(W*qr))';
    dFdee = -N/2*P.*P';
 
    % dF/dhh = -trace{P*Q{i}*P*Q{j}}
    %----------------------------------------------------------------------
    dhdh  = dedh(rel,as)*diag(exp(h(as)));
    dFdh  = dhdh'*dFde;
    dFdhh = dhdh'*dFdee*dhdh;
    
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP(as,as)*e(as);
    dFdhh = dFdhh - hP(as,as);
 
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = spm_dx(dFdhh,dFdh)/log(k + 2);
    h(as) = h(as) + dh;
 
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    % bar(h),drawnow
    dF      = dFdh'*dh;
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(dF));
    
    % and ARD
    %----------------------------------------------------------------------
    if dF < 1e-2 || k == K
        break
    else
        as             = h > -16;
        h(~as)         = hE(~as);
        h(find(h > 1)) = 1;
    end
 
end
 
 
% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
Ph  = hP;
Ph(as,as) = Ph(as,as) - dFdhh;
if nargout > 3
    
    % tr(hP*inv(Ph)) - nh + tr(pP*inv(Pp)) - np (pP = 0)
    %----------------------------------------------------------------------
    Ft = trace(hP*inv(Ph)) - length(Ph) - length(Cq);
            
    % complexity - KL(Ph,hP)
    %----------------------------------------------------------------------
    Fc = Ft/2 + e'*hP*e/2 + spm_logdet(Ph*inv(hP))/2 - N*spm_logdet(Cq)/2;
    
    % Accuracy - ln p(Y|h)
    %----------------------------------------------------------------------
    Fa = Ft/2 - trace(C*U*YY*U')/2 - N*n*log(2*pi)/2 - N*spm_logdet(C)/2;
    
    % Free-energy
    %----------------------------------------------------------------------
    F  = Fa - Fc - N*n*log(sY)/2;
    
end
 
% return exp(h) if log-normal hyperpriors
%--------------------------------------------------------------------------
h  = sY*exp(h)./sh;
C  = sY*C;

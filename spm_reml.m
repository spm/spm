function [V,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N);
% ReML estimation of [improper] covariance components from y*y'
% FORMAT [C,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of h
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates.
%
% see also: spm_reml_sc for the equivalent scheme using log-normal
% hyperpriors
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
 
% John Ashburner & Karl Friston
% $Id: spm_reml.m 3527 2009-11-02 20:27:13Z karl $
 
% assume a single sample if not specified
%--------------------------------------------------------------------------
try, N; catch, N  = 1;  end
 
% default number of iterations
%--------------------------------------------------------------------------
try, K; catch, K  = 128; end
 
% catch NaNs
%--------------------------------------------------------------------------
W     = Q;
q     = find(all(isfinite(YY)));
YY    = YY(q,q);
for i = 1:length(Q)
    Q{i} = Q{i}(q,q);
end
 
% initialise h
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);

% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
else
    X = spm_svd(X(q,:));
end

% initialise and specify hyperpriors
%==========================================================================
for i = 1:m
    h(i,1) = any(diag(Q{i}));
end
hE  = sparse(m,1);
hP  = speye(m,m)/exp(32);
dF  = Inf;
t   = 4;
 
 
% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:K
 
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        C = C + Q{i}*h(i);
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
 
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    P     = iC - iCX*Cq*iCX';
    U     = speye(n) - P*YY/N;
    for i = 1:m
 
        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        PQ{i}     = P*Q{i};
        dFdh(i,1) = -sum(sum(PQ{i}'.*U))*N/2;
 
    end
 
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m
 
            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j) = -sum(sum(PQ{i}'.*PQ{j}))*N/2;
            dFdhh(j,i) =  dFdhh(i,j);
 
        end
    end
    
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;
 
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = spm_dx(dFdhh,dFdh,{t});
    h     = h + dh;
 
    % Convergence (1% change in log-evidence)
    %======================================================================
    
    % update regulariser
    %----------------------------------------------------------------------
    if dF < dFdh'*dh, t = max(t - 1,-4); end
    dF    = dFdh'*dh;
    fprintf('%-24s: %3i %14s%e [%i]\n','  ReML Iteration',k,'...',full(dF),t);
    
    % final estimate of covariance (with missing data points)
    %----------------------------------------------------------------------
    if dF < 1e-1
        V     = 0;
        for i = 1:m
            V = V + W{i}*h(i);
        end
        break
    end
 
end
 
% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
Ph    = -dFdhh;
if nargout > 3
    
    % tr(hP*inv(Ph)) - nh + tr(pP*inv(Pp)) - np (pP = 0)
    %----------------------------------------------------------------------
    Ft = trace(hP*inv(Ph)) - length(Ph) - length(Cq);
    
    % complexity - KL(Ph,hP)
    %----------------------------------------------------------------------
    Fc = Ft/2 + e'*hP*e/2 + spm_logdet(Ph*inv(hP))/2 - N*spm_logdet(Cq)/2;
    
    % Accuracy - ln p(Y|h)
    %----------------------------------------------------------------------
    Fa = Ft/2 - trace(C*P*YY*P)/2 - N*n*log(2*pi)/2 - N*spm_logdet(C)/2;
    
    % Free-energy
    %----------------------------------------------------------------------
    F  = Fa - Fc;
    
end

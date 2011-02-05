function  [C,h,Ph,F,Fa,Fc] = spm_reml_A(YY,X,Q,N,hE,hC,V)
% ReML estimation of covariance components from y*y' - factored components
% FORMAT [C,h,Ph,F,Fa,Fc] = spm_reml_A(YY,X,Q,N,[hE,hC,V]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components (factors)
% N   - number of samples
%
% hE  - hyperprior expectation [default = 0]
% hC  - hyperprior covariance  [default = 256]
% V   - fixed covariance component
%
% C   - (m x m) estimated errors: C = A*A': A = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of h
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% Performs a Fisher-Scoring ascent on F to find MAP variance parameter
% estimates. NB: uses weakly informative normal hyperpriors on the factors.
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
% $Id: spm_reml_A.m 4195 2011-02-05 18:39:13Z ged $


% assume a single sample if not specified
%--------------------------------------------------------------------------
try, N; catch, N  = 1;      end
try, V; catch, V  = 1e-16;  end

% initialise h
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
h     = zeros(m,1) + exp(-4);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);

% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
else
    X = spm_svd(X,0);
end

% check fixed component
%--------------------------------------------------------------------------
if length(V) == 1
    V = V*speye(n,n);
end


% initialise and specify hyperpriors
%==========================================================================

% hyperpriors
%--------------------------------------------------------------------------
try, hE = hE(:);        catch, hE = 0;    end
try, hP = spm_inv(hC);  catch, hP = 1/256; end

% check sise
%--------------------------------------------------------------------------
if length(hE) < m, hE = hE(1)*ones(m,1);   end
if length(hP) < m, hP = hP(1)*speye(m,m);  end

% intialise h: so that sum(exp(h)) = 1
%--------------------------------------------------------------------------
if any(diag(hP) > exp(16))
    h = hE;
end

% ReML (EM/VB)
%--------------------------------------------------------------------------
dF    = Inf;
t     = 4;
for k = 1:32
    
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    A     = 0;
    for i = 1:m
        A = A + Q{i}*h(i);
    end
    C     = V + A*A';
    iC    = spm_inv(C);
    
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
    PQ = cell(m, 1);
    for i = 1:m
        
        % dF/dh
        %------------------------------------------------------------------
        PQ{i}   = P*(A*Q{i}' + Q{i}'*A);
        dFdh(i) = -sum(sum(PQ{i}.*U'))*N/2; % = -trace(PQ{i}*U)*N/2;
        
    end
    
    % Expected curvature E{dF/dhh} (second derivatives)
    % dF/dhh = -trace{P*Q{i}*P*Q{j}} NB trace(A*B)==sum(sum(A.*B'))
    %----------------------------------------------------------------------
    for i = 1:m
        dFdhh(i,i) = -sum(sum(PQ{i}.*PQ{i}'))*N/2;
        for j = i+1:m
            dFdhh(i,j) = -sum(sum(PQ{i}.*PQ{j}'))*N/2;
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
    
    
    % predicted change in F - increase regularisation if increasing
    %----------------------------------------------------------------------
    pF    = dFdh'*dh;
    if pF > dF
        t = t - 1;
    else
        t = t + 1/8;
    end
    dF    = pF;
    
    % convergence
    %----------------------------------------------------------------------
    fprintf('%s %-23d: %10s%e [%+3.2f]\n','  ReML Iteration',k,'...',full(dF),t);
    if dF < 1e-2
        break
    end
end

% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
Ph    = -dFdhh;
if nargout > 3
    
    % tr(hP*inv(Ph)) - nh + tr(pP*inv(Pp)) - np (pP = 0)
    %----------------------------------------------------------------------
    Ft = trace(hP/Ph) - length(Ph) - length(Cq);
    
    % complexity - KL(Ph,hP)
    %----------------------------------------------------------------------
    Fc = Ft/2 + e'*hP*e/2 + spm_logdet(Ph/hP)/2 - N*spm_logdet(Cq)/2;
    
    % Accuracy - ln p(Y|h)
    %----------------------------------------------------------------------
    Fa = Ft/2 - trace(C*P*YY*P)/2 - N*n*log(2*pi)/2  - N*spm_logdet(C)/2;
    
    % Free-energy
    %----------------------------------------------------------------------
    F  = Fa - Fc;
    
end

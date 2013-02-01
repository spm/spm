function [V,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N,D,t,hE,hP)
% ReML estimation of [improper] covariance components from y*y'
% FORMAT [C,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N,D,t,hE,hP);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
%
% N   - number of samples                 (default 1)
% D   - Flag for positive-definite scheme (default 0)
% t   - regularisation                    (default 4)
% hE  - hyperprior                        (default 0)
% hP  - hyperprecision                    (default 1e-16)
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
% $Id: spm_reml.m 5223 2013-02-01 11:56:05Z ged $
 
 
% check defaults
%--------------------------------------------------------------------------
try, N;  catch, N  = 1;     end       % assume a single sample if not specified
try, K;  catch, K  = 32;    end       % default number of iterations
try, D;  catch, D  = 0;     end       % default checking
try, t;  catch, t  = 4;     end       % default regularisation
try, hE; catch, hE = 0;     end       % default hyperprior
try, hP; catch, hP = 1e-16; end       % default hyperprecision
 
% catch NaNs
%--------------------------------------------------------------------------
W     = Q;
q     = find(all(isfinite(YY)));
YY    = YY(q,q);
for i = 1:length(Q)
    Q{i} = Q{i}(q,q);
end
 
% dimensions
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
 
% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
else
    X = spm_svd(X(q,:),0);
end
 
% initialise h and specify hyperpriors
%==========================================================================
h   = zeros(m,1);
for i = 1:m
    h(i,1) = any(diag(Q{i}));
end
hE  = sparse(m,1) + hE;
hP  = speye(m,m)*hP;
dF  = Inf;
D   = 8*(D > 0);
 
 
% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:K

    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        C = C + Q{i}*h(i);
    end
 
    % positive [semi]-definite check
    %----------------------------------------------------------------------
    for i = 1:D
        if min(real(eig(full(C)))) < 0

            % increase regularisation and re-evaluate C
            %--------------------------------------------------------------
            t     = t - 1;
            h     = h - dh;
            dh    = spm_dx(dFdhh,dFdh,{t});
            h     = h + dh;
            C     = sparse(n,n);
            for i = 1:m
                C = C + Q{i}*h(i);
            end
        else
            break
        end
    end


    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iC     = spm_inv(C);
    iCX    = iC*X;
    if ~isempty(X)
        Cq = spm_inv(X'*iCX);
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
        dFdh(i,1) = -spm_trace(PQ{i},U)*N/2;

    end

    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m

            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j) = -spm_trace(PQ{i},PQ{j})*N/2;
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
        t = t + 1/4;
    end
    
    % revert to SPD checking, if near phase-transition
    %----------------------------------------------------------------------
    if ~isfinite(pF) || abs(pF) > 1e6
        [V,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N,1,t - 2);
        return
    else
        dF = pF;
    end
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    fprintf('%s %-23d: %10s%e [%+3.2f]\n','  ReML Iteration',k,'...',full(pF),t);
 
    % final estimate of covariance (with missing data points)
    %----------------------------------------------------------------------
    if dF < 1e-1, break, end
end

 
% re-build predicted covariance
%==========================================================================
V     = 0;
for i = 1:m
    V = V + W{i}*h(i);
end
 
% check V is positive semi-definite (if not already checked)
%==========================================================================
if ~D
    if min(eig(V)) < 0
        [V,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N,1,2,hE(1),hP(1));
        return
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

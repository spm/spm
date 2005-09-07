function [C,h,Ph,F] = spm_reml(YY,X,Q,N,OPT);
% ReML estimation of covariance components from y*y'
% FORMAT [C,h,Ph,F] = spm_reml(YY,X,Q,N,[OPT]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x n) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% OPT(1) = 1 : log-normal hyperparameterisation
% OPT(2) = 1 : invoke priors on h - default [0 0]
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
% $Id: spm_reml.m 222 2005-09-07 16:49:37Z karl $
 
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
    OPT = [0 0];
end
 
% ortho-normalise X
%--------------------------------------------------------------------------
X     = sparse(orth(full(X)));
[n p] = size(X);
 
% ensure estimable components
%--------------------------------------------------------------------------
m     = length(Q);
dh    = sparse(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
for i = 1:m
    for j = i:m
        dFdhh(i,j) = sum(sum(Q{i}.*Q{j}'));
        dFdhh(j,i) = dFdhh(i,j);
    end
end
 
% enforce hyperpriors if curvature is rank deficient
%--------------------------------------------------------------------------
if rank(dFdhh) < m
    OPT(2) = 1;
end
 
if OPT(1)
 
    % log-normal hyperpriors - expectation and precision
    %----------------------------------------------------------------------
    hE    = ones(m,1);
    hP    = speye(m,m)/8;
    h     = hE;
    
else
 
    % Gaussian hyperpriors - expectation and precision
    %----------------------------------------------------------------------
    hE    = sparse(m,1);
    hP    = speye(m,m)/1e6;
 
    % initialise h
    %----------------------------------------------------------------------
    C     = [];
    R     = speye(n,n) - X*X';
    for i = 1:m
        RQR = R*Q{i}*R';
        C   = [C RQR(:)];
    end
    R     = R*YY*R'/N;
    h     = inv(C'*C)*(C'*R(:));
 
end
 
% hyperpriors
%--------------------------------------------------------------------------
if ~OPT(2)
    hP = hP*0;
end
 
% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:32
 
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        if OPT(1)
            C = C + Q{i}*exp(h(i));
        else
            C = C + Q{i}*h(i);
        end
    end
    iC    = inv(C);
 
    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iCX   = iC*X;
    Cq    = inv(X'*iCX);
    XCXiC = X*Cq*iCX';
     
    % M-step: ReML estimate of hyperparameters
    %======================================================================
 
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    P     = iC - iC*XCXiC;
    PYY   = P*YY;
    for i = 1:m
 
        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        PQ{i}     = P*Q{i};
        if OPT(1)
            PQ{i} = PQ{i}*exp(h(i));
        end
        dFdh(i)   = -trace(PQ{i})*N/2 + sum(sum(PQ{i}.*PYY'))/2;

 
    end
    dFdh  = dFdh - hP*(h - hE);
 
    % Expected curvature E{ddF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m
 
            % ddF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j)  = -sum(sum(PQ{i}.*PQ{j}'))*N/2 - hP(i,j);
            dFdhh(j,i)  =  dFdhh(i,j);

        end
    end
 
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = -pinv(dFdhh)*dFdh;
    h     = h + dh;
 
    % Convergence (1% change in log-evidence)
    %======================================================================
    w     = dFdh'*dh;
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(w));
    if w < 1e-2, break, end
 
end
 
% hyperpriors - precision
%--------------------------------------------------------------------------
if OPT(1)
    h  = exp(h);
end

% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
if nargout > 3

    % compute condotional covariance of h
    %----------------------------------------------------------------------
    for i = 1:m
        CP{i} = -Q{i}*iC;
        if OPT(1)
            CP{i} = CP{i}*exp(h(i));
        end
    end
    
    % P(h) = -ddF/dhh - ...
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m
            CPCP     =  CP{i}*CP{j};
            Ph(i,j)  = -sum(sum(CPCP.*XCXiC'))*N;
            Ph(j,i)  =  Ph(i,j);
        end
    end
    Ph = Ph - dFdhh;

    % log evidence = F
    %----------------------------------------------------------------------
    F = - trace(C*PYY*P)/2 ...
        - N*n*log(2*pi)/2 ...
        - N*spm_logdet(C)/2 ...
        + N*spm_logdet(Cq)/2 ...
        +   spm_logdet(hP)/2 ...
        -   spm_logdet(Ph)/2 ...
        - N*p/2 - m/2;
end





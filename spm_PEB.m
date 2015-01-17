function [C,P,F] = spm_PEB(y,P,OPT)
% parametric empirical Bayes (PEB) for hierarchical linear models
% FORMAT [C,P,F] = spm_PEB(y,P,OPT)
%
% y       - (n x 1)     response variable
%
% MODEL SPECIFICATION
%
% P{i}.X  - (n x m)     ith level design matrix i.e: constraints on <Eb{i - 1}>
% P{i}.C  - {q}(n x n)  ith level contraints on Cov{e{i}} = Cov{b{i - 1}}
%
% OPT    - enforces positively constraints on the covariance hyperparameters
%          by adopting a log-normal [flat] hyperprior. default = 0
%
% POSTERIOR OR CONDITIONAL ESTIMATES
%
% C{i}.E  - (n x 1)     conditional expectation E{b{i - 1}|y}
% C{i}.C  - (n x n)     conditional covariance  Cov{b{i - 1}|y} = Cov{e{i}|y}
% C{i}.M  - (n x n)     ML estimate of Cov{b{i - 1}} = Cov{e{i}}
% C{i}.h  - (q x 1)     ith level ReML  hyperparameters for covariance:
%                       Cov{e{i}} = P{i}.h(1)*P{i}.C{1} +  ...
%
% LOG EVIDENCE
%
% F       - [-ve] free energy F = log evidence = p(y|X,C)
%
% If P{i}.C is not a cell the covariance at that level is assumed to be kown
% and Cov{e{i}} = P{i}.C (i.e. the hyperparameter is fixed at 1)
%
% If P{n}.C is not a cell this is taken to indicate that a full Bayesian
% estimate is required where P{n}.X is the prior expectation and P{n}.C is
% the known prior covariance.  For consistency, with PEB, this is implemented
% by setting b{n} = 1 through appropriate constraints at level {n + 1}.
%
% To implement non-hierarchical Bayes with priors on the parameters use
% a two level model setting the second level design matrix to zeros.
%__________________________________________________________________________
%
% Returns the moments of the posterior p.d.f. of the parameters of a
% hierarchical linear observation model under Gaussian assumptions
%
%                            y = X{1}*b{1} + e{1}
%                         b{1} = X{2}*b{2} + e{2}
%                                 ...
%
%                     b{n - 1} = X{n}*b{n} + e{n}
%
% e{n} ~ N{0,Ce{n}}
%
% using Parametic Emprical Bayes (PEB)
%
% Ref: Dempster A.P., Rubin D.B. and Tsutakawa R.K. (1981) Estimation in
% covariance component models.  J. Am. Stat. Assoc. 76;341-353
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_PEB.m 6305 2015-01-17 12:40:51Z karl $

% set default
%--------------------------------------------------------------------------
try
    OPT;
catch
    OPT = 0;
end

% number of levels (p)
%--------------------------------------------------------------------------
M     = 32;                                  % maximum number of iterations
p     = length(P);

% check covariance constraints - assume i.i.d. errors conforming to X{i}
%--------------------------------------------------------------------------
for i = 1:p
    if ~isfield(P{i},'C')
        [n,m] = size(P{i}.X);
        if i == 1
            P{i}.C  = {speye(n,n)};
        else
            for j = 1:m
                k         = find(P{i}.X(:,j));
                P{i}.C{j} = sparse(k,k,1,n,n);
            end
        end
    end

end

% Construct augmented non-hierarchical model
%==========================================================================

% design matrix and indices
%--------------------------------------------------------------------------
I     = {0};
J     = {0};
K     = {0};
XX    = [];
X     = 1;
for i = 1:p

    % design matrix
    %----------------------------------------------------------------------
    X     = X*P{i}.X;
    XX    = [XX X];

    % indices for ith level parameters
    %----------------------------------------------------------------------
    [n,m] = size(P{i}.X);
    I{i}  = (1:n) + I{end}(end);
    J{i}  = (1:m) + J{end}(end);

end

% augment design matrix and data
%--------------------------------------------------------------------------
n        = size(XX,2);
XX       = [XX; speye(n,n)];
y        = [y; sparse(n,1)];

% last level constraints
%--------------------------------------------------------------------------
n        = size(P{p}.X,2);
I{p + 1} = (1:n) + I{end}(end);
q        = I{end}(end);
Cb       = sparse(q,q);
if ~iscell(P{end}.C)

    % Full Bayes: (i.e. Cov(b) = 0, <b> = 1)
    %----------------------------------------------------------------------
    y( I{end})        = sparse(1:n,1,1);
else

    % Empirical Bayes: uniform priors (i.e. Cov(b) = Inf, <b> = 0)
    %----------------------------------------------------------------------
    Cb(I{end},I{end}) = sparse(1:n,1:n,exp(32));
end


% assemble augmented  constraints Q: Cov{e} = Cb + h(i)*Q{i} + ...
%==========================================================================
if ~isfield(P{1},'Q')

    % covariance contraints Q on Cov{e{i}} = Cov{b{i - 1}}
    %----------------------------------------------------------------------
    h     = [];
    Q     = {};
    for i = 1:p

        % collect constraints on prior covariances - Cov{e{i}}
        %------------------------------------------------------------------
        if iscell(P{i}.C)
            m     = length(P{i}.C);
            for j = 1:m
                [u,v,s]    = find(P{i}.C{j});
                u          = u + I{i}(1) - 1;
                v          = v + I{i}(1) - 1;
                Q{end + 1} = sparse(u,v,s,q,q);
            end

            % indices for ith-level hyperparameters
            %--------------------------------------------------------------
            try
                K{i}  = (1:m) + K{end}(end);
            catch
                K{i}  = (1:m);
            end

        else

            % unless they are known - augment Cb
            %--------------------------------------------------------------
            [u,v,s] = find(P{i}.C + speye(length(P{i}.C))*1e-6);
            u       = u + I{i}(1) - 1;
            v       = v + I{i}(1) - 1;
            Cb      = Cb + sparse(u,v,s,q,q);

            % indices for ith-level hyperparameters
            %--------------------------------------------------------------
            K{i}  = [];

        end

    end

    % note overlapping bases - requiring 2nd order M-Step derivatives
    %----------------------------------------------------------------------
    m     = length(Q);
    d     = sparse(m,m);
    for i = 1:m
        XQX{i} = XX'*Q{i}*XX;
    end
    for i = 1:m
        for j = i:m
            o      = nnz(XQX{i}*XQX{j});
            d(i,j) = o;
            d(j,i) = o;
        end
    end

    % log-transform and save
    %----------------------------------------------------------------------
    h   = zeros(m,1);
    if OPT
        hP  = speye(m,m)/16;
    else
        hP  = speye(m,m)/exp(16);
        for i = 1:m
            h(i) = any(diag(Q{i}));
        end
    end
    P{1}.hP = hP;
    P{1}.Cb = Cb;
    P{1}.Q  = Q;
    P{1}.h  = h;
    P{1}.K  = K;
    P{1}.d  = d;

end
hP    = P{1}.hP;
Cb    = P{1}.Cb;
Q     = P{1}.Q;
h     = P{1}.h;
K     = P{1}.K;
d     = P{1}.d;

% Iterative EM
%--------------------------------------------------------------------------
m     = length(Q);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
for k = 1:M

    % inv(Cov(e)) - iC(h)
    %----------------------------------------------------------------------
    Ce    = Cb;
    for i = 1:m
        if OPT
            Ce = Ce + Q{i}*exp(h(i));
        else
            Ce = Ce + Q{i}*h(i);
        end
    end
    iC    = spm_inv(Ce,exp(-16));

    % E-step: conditional mean E{B|y} and covariance cov(B|y)
    %======================================================================
    iCX   = iC*XX;
    Cby   = spm_inv(XX'*iCX);
    B     = Cby*(iCX'*y);


    % M-step: ReML estimate of hyperparameters (if m > 0)
    %======================================================================
    if m == 0, break, end

    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    Py    = iC*(y - XX*B);
    iCXC  = iCX*Cby;
    for i = 1:m

        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        PQ{i}     = iC*Q{i} - iCXC*(iCX'*Q{i});
        if OPT
            PQ{i} = PQ{i}*exp(h(i));
        end
        dFdh(i)   = -trace(PQ{i})/2 + y'*PQ{i}*Py/2;

    end

    % Expected curvature E{ddF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m
            if d(i,j)

                % ddF/dhh = -trace{P*Q{i}*P*Q{j}}
                %----------------------------------------------------------
                dFdhh(i,j)  = -spm_trace(PQ{i},PQ{j})/2;
                dFdhh(j,i)  =  dFdhh(i,j);

            end
        end
    end

    % add hyperpriors
    %----------------------------------------------------------------------
    dFdhh = dFdhh - hP;
    
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = -pinv(dFdhh)*dFdh;
    h     = h + dh;

    % Convergence
    %======================================================================
    w     = norm(dh,1);
    
    fprintf('%-30s: %i %30s%e\n','  PEB Iteration',k,'...',full(w));
    
    % if dF < 0.01
    %----------------------------------------------------------------------
    if dFdh'*dh < 1e-2,     break, end
    
    % if dh^2 < 1e-8
    %----------------------------------------------------------------------
    if w < 1e-4,            break, end
    
    % if log-normal hyperpriors and h < exp(-16)
    %----------------------------------------------------------------------
    if OPT && all(h < -16), break, end

end

% place hyperparameters in P{1} and output structure for {n + 1}
%--------------------------------------------------------------------------
P{1}.h         = h + exp(-32);
C{p + 1}.E     = B(J{p});
C{p + 1}.M     = Cb(I{end},I{end});

% recursive computation of conditional means E{b|y}
%--------------------------------------------------------------------------
for i = p:-1:2
    C{i}.E     = B(J{i - 1}) + P{i}.X*C{i + 1}.E;
end

% hyperpriors - precision
%--------------------------------------------------------------------------
if OPT
    h          = exp(h);
end

% conditional covariances Cov{b|y} and ReML esimtates of Ce{i) = Cb{i - 1}
%--------------------------------------------------------------------------
for i = 1:p
    C{i + 1}.C = Cby(J{i},J{i});
    C{i}.M     = Ce(I{i},I{i});
    C{i}.h     = h(K{i});
end

% log evidence = ln p(y|X,C) = F = [-ve] free energy
%--------------------------------------------------------------------------
if nargout > 2

    % condotional covariance of h
    %----------------------------------------------------------------------
    Ph = -dFdhh;

    % log evidence = F
    %----------------------------------------------------------------------
    F = - Py'*Ce*Py/2 ...
        - length(I{1})*log(2*pi)/2 ...
        - spm_logdet(Ce)/2 ...
        - spm_logdet(Ph)/2 ...
        + spm_logdet(hP)/2 ...
        + spm_logdet(Cby)/2;
end

% warning
%--------------------------------------------------------------------------
if k == M, warning('maximum number of iterations exceeded'), end


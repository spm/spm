function [p0,X,F,f,NESS] = spm_ness_GN(M)
% nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT [p0,X,F,f,NESS] = spm_ness_GN(M)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.f   - dx/dt   = f(x,u,P)  {function string or m-file}
%    M.pE  - P       = parameters of equation of motion
%    M.x   - (n x 1) = x(0) = expansion point
%    M.W   - (n x n) - precision matrix of random fluctuations
% x    - cell array of vectors specifying evaluation grid
%
% p0      - nonequilibrium steady-state
% X       - evaluation points of state space
% F       - expected flow
% f       - original flow
%
% NESS.H  - expected Hessian
% NESS.J  - expected Jacobian
% NESS.E  - Lyapunov exponents
% NESS.H2 - expected Euclidean norm of Hessian
% NESS.J2 - expected Euclidean norm of Jacobian
% NESS.D2 - correlation dimension
% NESS.bS - p0 = spm_softmax(spm_dctmtx(nx,nb)*bS);
% NESS.nb - number of basis functions
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_hd.m 8000 2020-11-03 19:04:17QDb karl $

% scale space search
%==========================================================================
for s = 1:1
    
    % affine scaling of sample points
    %----------------------------------------------------------------------
    x   = mean(M.X);
    M.X = bsxfun(@minus,M.X,x)/s;
    M.X = bsxfun(@plus,M.X,x);
    
    % get domain of phase-space and polynomial basis set
    %======================================================================
    U      = spm_ness_U(M);       % basis set for this domain
    [~,~,~,o] = spm_polymtx(U.x,U.K);
    [nX,n] = size(U.X);
    

    if s == 1
        
        % initialise priors using iterated maximum likelihood solution
        %------------------------------------------------------------------
        [~,X,~,f,NESS] = spm_ness_hd(M);
        
        pE.Sp = NESS.bS;
        pE.Qp = NESS.bQ;
        pC.Sp = spm_zeros(pE.Sp) + exp(8);
        pC.Qp = spm_zeros(pE.Qp) + exp(8);
        pC    = diag(spm_vec(pC));
        
    end
    
    % constraints on solenoidal coefficients
    %----------------------------------------------------------------------
    r    = o < 3;
    for i = 1:n
        for j = (i + 1):n
            r = [r (o < 2)];
        end
    end
    R    = diag(~r);
    pC   = R*pC*R;

    % rotate priors to orthonormal coefficients
    %----------------------------------------------------------------------
    nu    = (n^2 - n)/2 + 1;
    u     = kron(eye(nu,nu),inv(U.v));
    pE    = spm_unvec(u*spm_vec(pE),pE);
    pC    = u*pC*u';


    % model specification
    %======================================================================
    S.Nmax = 32;                    % maximum number of iterations
    S.G    = @spm_NESS_gen;        % generative function
    S.FS   = @(y)y(:);             % generative function
    S.pE   = pE;                   % prior expectations (parameters)
    S.pC   = pC;                   % prior covariances  (parameters)
    S.hE   = 4;                    % prior expectation  (log-precision)
    S.hC   = 0;                    % prior covariances  (log-precision)
    
    % model inversion with Variational Laplace (Gauss Newton)
    %======================================================================
    [Ep,Cp] = spm_nlsi_GN(S,U,U.f');
    
    % emprical priors (Bayesian belief updating)
    %----------------------------------------------------------------------
    v   = kron(eye(nu,nu),U.v);
    pE  = spm_unvec(v*spm_vec(Ep),Ep);
    pC  = v*Cp*v';

end


% Bayesian model reduction
%--------------------------------------------------------------------------
% DCM.M.pE  = v*spm_vec(S.pE);
% DCM.M.pC  = v*S.pC*v;
% DCM.Ep    = v*spm_vec(Ep);
% DCM.Cp    = v*Cp*v;
% 
% RCM = spm_dcm_bmr_all(DCM,'All','BMS')

% recompute Hessian D*D*S
%--------------------------------------------------------------------------
H     = zeros(n,n,nX);
HH    = cell(n,n);
for i = 1:n
    for j = 1:n
        HH{i,j} = U.D{i}*U.b'*U.D{j}*Ep.Sp;
    end
end
for i = 1:n
    for j = 1:n
        for k = 1:nX
            H(i,j,k) = HH{i,j}(k);
        end
    end
end

% update NESS structure
%--------------------------------------------------------------------------
p0      = spm_softmax(U.b*Ep.Sp);         % NESS density
NESS.H  = spm_dot(H   ,p0);               % expected Hessian
NESS.H2 = spm_dot(H.^2,p0);               % expected Euclidean norm of Hessian
NESS.Ep = pE;                             % posterior expectations
NESS.Cp = pC;                             % posterior covariances
NESS.bS = pE.Sp;                          % p0 = spm_softmax(spm_polymtx(x,nb)*bS);
NESS.bQ = pE.Qp;                          % parameters of solenoidal operator


% reshape nonequilibrium steady-state density
%--------------------------------------------------------------------------
p0  = reshape(p0,U.nx);

% predicted flow: F   = Q*Dx*S - L
%--------------------------------------------------------------------------
F   = spm_NESS_gen(Ep,M,U);

return

function [NESS] = spm_ness_Lap(M,x)
% Nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT [NESS] = spm_ness_Lap(M,x)
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

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% get domain of phase-space and polynomial basis set
%==========================================================================
U      = spm_ness_U(M,x);      % get state space and flow
f      = U.f';                 % target flow
[nX,n] = size(U.X);            % number of sample points and states 

% get constraints (allowing G to be free parameters)
%--------------------------------------------------------------------------
J      = any(U.J,3);
J      = J'.*J;
[ks,kq,kg] = spm_NESS_constraints(U.o,J,M.K,M.L);
k      = find(~(kq | kg));
k      = find(~(kg));
nb     = size(U.b,2);
nQ     = n*n/2 + n/2;

pE.Qp  = zeros(nb*nQ,1);       % polynomial coefficients for solenoidal operator
pE.Rp  = zeros(1 ,n);          % polynomial coefficients for mean
pE.Sp  = zeros(nb,n,n);        % polynomial coefficients for Kernel

pC.Qp  = pE.Qp;
pC.Rp  = pE.Rp;
pC.Sp  = pE.Sp;

% constraints on solenoidal operator
%--------------------------------------------------------------------------
pC.Qp(k) = 1;

% constraints on mean
%--------------------------------------------------------------------------
for i = 1:n
    pE.Rp(1,i)  = mean(x{i});
    pC.Rp(1,i)  = 1;
end

% constraints on potential
%--------------------------------------------------------------------------
k     = sum(U.o) < M.K;
for i = 1:n
    for j = 1:n
        if i == j
           pE.Sp(1,i,j) = 1/8;
        end
        if j >= i
           pC.Sp(k,i,j) = 1;
        end
    end
end

% priors
%--------------------------------------------------------------------------
pC     = diag(spm_vec(pC));

% model specification
%==========================================================================
S.Nmax = 64;                   % maximum number of iterations
S.G    = @spm_NESS_gen_lap;    % generative function
S.FS   = @(y)y(:);             % generative function
S.pE   = pE;                   % prior expectations (parameters)
S.pC   = pC;                   % prior covariances  (parameters)
S.hE   = 8;                    % prior expectation  (log-precision)
S.hC   = 1/512;                % prior covariances  (log-precision)

% model inversion with Variational Laplace (Gauss Newton)
%==========================================================================
Ep    = spm_nlsi_GN(S,U,f);

% NESS density and expected flow
%--------------------------------------------------------------------------
[F,S,Q,L,H] = spm_NESS_gen_lap(Ep,M,U);           
p0          = spm_softmax(-S(:));

% update NESS structure
%--------------------------------------------------------------------------
NESS.H  = spm_dot(H   ,p0);           % expected Hessian
NESS.H2 = spm_dot(H.^2,p0);           % expected Euclidean norm of Hessian
NESS.pE = pE;                         % posterior expectations
NESS.Ep = Ep;                         % posterior expectations

NESS.p0 = reshape(p0,U.nx);           % nonequilibrium steady-state
NESS.X  = U.X;                        % evaluation points of state space
NESS.F  = F;                          % expected flow
NESS.f  = f;                          % original flow

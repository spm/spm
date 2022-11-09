function [NESS] = spm_ness_GN(M,x)
% Nonequilibrium steady-state under a Helmholtz decomposition
% FORMAT [NESS] = spm_ness_GN(M,x)
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

% initialise priors using iterated maximum likelihood solution
%--------------------------------------------------------------------------
NESS   = spm_ness_hd(M,x);     % estimate and fixed G
pE     = NESS.Ep;

% get constraints (allowing G to be free parameters)
%--------------------------------------------------------------------------
[ks,kq,kg] = spm_NESS_constraints(U.o,any(U.J,3),M.K,M.L);
k      = [ks (kq | kg)];
i      = find(~k);
np     = numel(k);
pC     = sparse(i,i,16,np,np);

% model specification
%==========================================================================
S.Nmax = 32;                   % maximum number of iterations
S.G    = @spm_NESS_gen;        % generative function
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
p0    = spm_softmax(U.b*Ep.Sp);
F     = spm_NESS_gen(Ep,M,U);           

% Hessian D*D*S
%--------------------------------------------------------------------------
H     = zeros(n,n,nX);
HH    = cell(n,n);
for i = 1:n
    for j = 1:n
        HH{i,j} = U.H{i,j}*Ep.Sp;
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
NESS.H  = spm_dot(H   ,p0);           % expected Hessian
NESS.H2 = spm_dot(H.^2,p0);           % expected Euclidean norm of Hessian
NESS.pE = pE;                         % posterior expectations
NESS.Ep = Ep;                         % posterior expectations

NESS.p0 = reshape(p0,U.nx);           % nonequilibrium steady-state
NESS.X  = U.X;                        % evaluation points of state space
NESS.F  = F;                          % expected flow
NESS.f  = f;                          % original flow

function [F] = spm_NESS_F(P,M)
% Generate flow (f) at locations (U.X)
% FORMAT [F,S,Q,L,H,D] = spm_NESS_gen(P,M)
% FORMAT [F,S,Q,L,H,D] = spm_NESS_gen(P,M,U)
% FORMAT [F,S,Q,L,H,D] = spm_NESS_gen(P,M,X)
%--------------------------------------------------------------------------
% P.Qp    - polynomial coefficients for solenoidal operator
% P.Sp    - polynomial coefficients for potential
%
% F       - polynomial approximation to flow
% S       - negative potential (log NESS density)
% Q       - flow operator (R + G) with solenoidal and symmetric parts
% L       - correction term for derivatives of solenoidal flow
% H       - Hessian
% D       - potential gradients
%
% U = spm_ness_U(M)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.X  - sample points
%    M.W  - (n x n) - precision matrix of random fluctuations
%    M.K  - order of polynomial expansion
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


% get basis or expansion from M.X (or M.x)
%--------------------------------------------------------------------------

% get basis set and derivatives
%----------------------------------------------------------------------
U    = spm_ness_U(M);
    
% dimensions and correction terms to flow operator
%==========================================================================
n     = numel(U.D);

% predicted flow: F   = Q*D*S - L
%--------------------------------------------------------------------------
DS    = zeros(n,1);
for j = 1:n
    DS(j) = U.D{j}*P.Sp;
end

% quadratic forms
%==========================================================================
ih = [2,3];
is = 1;
ia = 4;
im = [5,6];
ip = [1,4,5,6];                         % particular states
in = [2:n];
x  = M.X';

% generative model (joint density)
%--------------------------------------------------------------------------
[m,C]  = spm_ness_cond(n,3,P.Sp);
Pi     = inv(C);
Ji     = (x - m)'*Pi*(x - m)/2;
dJidx  = Pi*(x - m)

% marginal over particular states (free energy)
%--------------------------------------------------------------------------
Pj     = inv(C(ip,ip));
Jj     = (x(ip) - m(ip))'*Pj*(x(ip) - m(ip))/2;
dJjdp  = Pj*(x(ip) - m(ip))

% explicit form for free energy
%==========================================================================

% likelihood
%--------------------------------------------------------------------------
[ml,Cl] = spm_ness_cond(n,3,P.Sp,in,x(in));
Pl      = inv(Cl);

% Prior
%--------------------------------------------------------------------------
mp      = m(in);
Pp      = inv(C(in,in));

F       = (x(in) - mp)'*Pp*(x(in) - mp)/2;
dFdp    = Pp*(x(in) - mp)

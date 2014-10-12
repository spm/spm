function [s,u] = spm_ssm2s(P,M,TOL)
% Converts state-space (M) representation to eigenspectrum
% FORMAT [s,u] = spm_ssm2s(P,M)
%
% P    - model parameters
% M    - model (with flow M.f and expansion point M.x and M.u)
% TOL  - optional upper bound for  principality exponent  (default -4)
%
% S    - (sorted) eigenspectrum or Lyapunov exponents
% V    - associated eigenvectors
%
% csd  - cross spectral density
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ssm2s.m 6233 2014-10-12 09:43:50Z karl $

% preliminaries
%--------------------------------------------------------------------------
if nargin < 3, TOL = -4; end


% Steady state solution
%--------------------------------------------------------------------------
M.x = spm_dcm_neural_x(P,M);

try
    M.u;
catch
    M.u = zeros(M.l,1);
end

% Jacobian and delay operator - if not specified already
%--------------------------------------------------------------------------
if nargout(M.f) >= 3
    [f,dfdx,D] = feval(M.f,M.x,M.u,P,M);
    
elseif nargout(M.f) == 2
    [f,dfdx]   = feval(M.f,M.x,M.u,P,M);
    D          = 1;
else
    dfdx       = spm_diff(M.f,M.x,M.u,P,M,1);
    D          = 1;
end

dfdx   = D*dfdx;
dfdu   = D*spm_diff(M.f,M.x,M.u,P,M,2);
[u,s]  = eig(full(dfdx),'nobalance');
s      = diag(s);

% eigenvectors
%--------------------------------------------------------------------------
u      = u*diag(pinv(u)*dfdu);

% condition slow eigenmodes
%--------------------------------------------------------------------------
s      = 1j*imag(s) + min(real(s),TOL);

% principal eigenmodes (highest imaginary value)
%--------------------------------------------------------------------------
[d,i]  = sort(imag(s),'descend');
u      = u(:,i);
s      = s(i);

% principal eigenmodes (highest real value)
%--------------------------------------------------------------------------
j      = find(~imag(s));
[d,i]  = sort(real(s(j)),'descend');
u(:,j) = u(:,j(i));
s(j)   = s(j(i));




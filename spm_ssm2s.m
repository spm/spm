function [S] = spm_ssm2s(P,M,U)
% Converts state-space (M) representation to eigenspectrum
% FORMAT [S] = spm_ssm2s(P,M,U)
%
% P    - model parameters
% M    - model (with flow M.f and expansion point M.x and M.u)
% U    - induces expansion about steady state
%
% S    - (SPM) eigenspectrum
% Hz   - vector of frequencies (Hz)
%
% csd  - cross spectral density
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ssm2s.m 5448 2013-04-25 11:08:52Z guillaume $


% Steady state solution
%--------------------------------------------------------------------------
try
    J  = M.J;
catch
    J  = [2; 4];
end


% Steady state solution
%--------------------------------------------------------------------------
if nargin > 2
    M.x = spm_dcm_neural_x(P,M);
end

% Jacobian and delay operator - if not specified already
%--------------------------------------------------------------------------
if nargout(M.f) == 3
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
[u,s]  = eig(full(dfdx));
s      = diag(s);

% alert to unstable eigenmodes
%--------------------------------------------------------------------------
if any(real(s) > 0)
    disp('warning unstable modes')
end

% remove unstable eigenmodes
%--------------------------------------------------------------------------
i     = find(real(s) < 0);
s     = s(i);
u     = u(:,i);

% condition slow eigenmodes
%--------------------------------------------------------------------------
s     = 1j*imag(s) - log(32 + exp(real(-s)));


% number of imaginary eigenmodes
%--------------------------------------------------------------------------
try
    m = M.Nm;
catch
    m = sum(imag(s) > 0);
end

% principal eigenmodes (highest imaginary value)
%--------------------------------------------------------------------------
[q,i]  = sort(imag(s),'descend');
u      = u(:,i);
s      = s(i);

% principal eigenmodes (highest real value)
%--------------------------------------------------------------------------
j      = find(~imag(s));

[r,i]  = sort(real(s(j)),'descend');
u(:,j) = u(:,j(i));
s(j)   = s(j(i));

% weights
%--------------------------------------------------------------------------
w      = u*diag(pinv(u)*dfdu);
w      = abs(w(J,1:m));
w      = w*32/max(w(:));


% principal eigenmodes (most unstable)
%--------------------------------------------------------------------------
S      =  [s(1:m); spm_vec(w)];

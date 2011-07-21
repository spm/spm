function [M0,M1,L1,L2] = spm_bireduce(M,P,Q)
% reduction of a fully nonlinear MIMO system to Bilinear form
% FORMAT [M0,M1,L1,L2] = spm_bireduce(M,P,D);
%
% M   - model specification structure
% Required fields:
%   M.f   - dx/dt    = f(x,u,P,M)                 {function string or m-file}
%   M.g   - y(t)     = g(x,u,P,M)                 {function string or m-file}
%   M.bi  - bilinear form [M0,M1,L1,L2] = bi(M,P) {function string or m-file}
%   M.m   - m inputs
%   M.n   - n states
%   M.l   - l outputs
%   M.x   - (n x 1) = x(0) = expansion point: defaults to x = 0;
%   M.u   - (m x 1) = u    = expansion point: defaults to u = 0;
%
% P   - model parameters
% D   - delay operator df/dx -> D*df/dx (default D = 1)
%
% A Bilinear approximation is returned where the states are
%
%        q(t) = [1; x(t) - x(0)]
%
%___________________________________________________________________________
% Returns Matrix operators for the Bilinear approximation to the MIMO
% system described by
%
%       dx/dt = f(x,u,P)
%        y(t) = g(x,u,P)
%
% evaluated at x(0) = x and u = 0
%
%       dq/dt = M0*q + u(1)*M1{1}*q + u(2)*M1{2}*q + ....
%        y(i) = L1(i,:)*q + q'*L2{i}*q/2;
%
%--------------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_bireduce.m 4401 2011-07-21 12:34:46Z karl $


% set up
%==========================================================================
try
    Q;
catch
    Q = 1;
end

% create inline functions
%--------------------------------------------------------------------------
try
    funx = fcnchk(M.f,'x','u','P','M');
catch
    M.f  = inline('sparse(0,1)','x','u','P','M');
    M.n  = 0;
    M.x  = sparse(0,0);
    funx = fcnchk(M.f,'x','u','P','M');
end

% expansion pointt
%--------------------------------------------------------------------------
x     = M.x;           
try
    u = spm_vec(M.u);
catch
    u = sparse(M.m,1);
end

% Partial derivatives for 1st order Bilinear operators
%==========================================================================

% f(x(0),0) and derivatives
%--------------------------------------------------------------------------
[dfdxu dfdx] = spm_diff(funx,x,u,P,M,[1 2]);
[dfdu  f0]   = spm_diff(funx,x,u,P,M,2);
f0           = spm_vec(f0);
m            = length(dfdxu);          % m inputs
n            = length(f0);             % n states

% delay operator
%--------------------------------------------------------------------------
f0    = Q*f0;
dfdx  = Q*dfdx;
dfdu  = Q*dfdu;
for i = 1:m
    dfdxu{i} = Q*dfdxu{i};
end


% Bilinear operators
%==========================================================================

% Bilinear operator - M0
%--------------------------------------------------------------------------
M0    = spm_cat({0                     []    ;
                (f0 - dfdx*spm_vec(x)) dfdx});

% Bilinear operator - M1 = dM0/du
%--------------------------------------------------------------------------
M1    = {};
for i = 1:m
    M1{i} = spm_cat({0,                                []        ;
                    (dfdu(:,i) - dfdxu{i}*spm_vec(x)), dfdxu{i}});
end

if nargout < 3, return, end


% Output operators
%==========================================================================

% add observer if not specified
%--------------------------------------------------------------------------
try
    fung = fcnchk(M.g,'x','u','P','M');
catch
    M.g  = inline('spm_vec(x)','x','u','P','M');
    M.l  = n;
    fung = fcnchk(M.g,'x','u','P','M');
end

% g(x(0),0)
%--------------------------------------------------------------------------
[dgdx g0] = spm_diff(fung,x,u,P,M,1);
g0        = spm_vec(g0);
l         = length(g0);

% Output matrices - L1
%--------------------------------------------------------------------------
L1    = spm_cat({(spm_vec(g0) - dgdx*spm_vec(x)), dgdx});


if nargout < 4, return, end

% Output matrices - L2
%--------------------------------------------------------------------------
dgdxx = spm_diff(fung,x,u,P,M,[1 1]);
for i = 1:l
    for j = 1:n
        D{i}(j,:) = dgdxx{j}(i,:);
    end
end
for i = 1:l
    L2{i} = spm_cat(spm_diag({0, D{i}}));
end
    

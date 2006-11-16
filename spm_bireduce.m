function [M0,M1,L1,L2] = spm_bireduce(M,P)
% reduction of a fully nonlinear MIMO system to Bilinear form
% FORMAT [M0,M1,L1,L2] = spm_bireduce(M,P);
%
% M   - model specification structure
% Required fields:
%   M.f   - dx/dt    = f(x,u,P)                     {function string or m-file}
%   M.g   - y(t)     = l(x,u,P)                     {function string or m-file}
%   M.bi  - bilinear form [M0,M1,L1,L2] = bi(M,P)   {function string or m-file}
%   M.m   - m inputs
%   M.n   - n states
%   M.l   - l outputs
%   M.x   - (n x 1) = x(0) = expansion point
%
% P   - model parameters
%
% A Bilinear approximation is returned where the states are
%
%		 q(t) = [1; x(t) - x(0)]
%
%___________________________________________________________________________
% Returns Matrix operators for the Bilinear approximation to the MIMO
% system described by
%
%		dx/dt = f(x,u,P)
% 		 y(t) = g(x,u,P)
%
% evaluated at x(0) = x and u = 0
%
%		dq/dt = M0*q + u(1)*M1{1}*q + u(2)*M1{2}*q + ....
%		 y(i) = L1(i,:)*q + q'*L2{i}*q/2;
%
%--------------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_bireduce $


% set up
%==========================================================================

% add [0] states if not specified
%--------------------------------------------------------------------------
if ~isfield(M,'f')
    M.f = inline('sparse(0,1)','x','u','P');
    M.n = 0;
    M.x = sparse(0,0);
end

% sizes and input
%--------------------------------------------------------------------------
m     = M.m;					% m inputs
n     = M.n;					% n states
l     = M.l;					% l ouputs
x     = M.x;                    % expansion point

try
    u = spm_vec(M.u);
catch
    u = sparse(m,1);
end

% create inline functions
%--------------------------------------------------------------------------
funx  = fcnchk(M.f,'x','u','P');
funl  = fcnchk(M.g,'x','u','P');


% f(x(0),0) and l(x(0),0)
%--------------------------------------------------------------------------
f0    = feval(funx,x,u,P);
l0    = feval(funl,x,u,P);


% Partial derivatives for 1st order Bilinear operators
%==========================================================================

% derivatives
%--------------------------------------------------------------------------
dfdx  = spm_diff(funx,x,u,P,1);
dfdxu = spm_diff(funx,x,u,P,[1 2]);
dfdu  = spm_diff(funx,x,u,P,2);

% Bilinear operators
%==========================================================================

% Bilinear operator - M0
%--------------------------------------------------------------------------
M0    = spm_cat({0                      []    ;
                 (f0 - dfdx*spm_vec(x)) dfdx});

% Bilinear operator - M1 = dM0/du
%--------------------------------------------------------------------------
for i = 1:m
    M1{i} = spm_cat({0,                               []         ;
                    (dfdu(:,i) - dfdxu{i}*spm_vec(x)), dfdxu{i}});
end

if nargout < 3, return, end

% Output matrices - L1
%--------------------------------------------------------------------------
dldx  = spm_diff(funl,x,u,P,1);
L1    = spm_cat({(l0 - dldx*spm_vec(x)), dldx});

if nargout < 4, return, end

% Output matrices - L2
%--------------------------------------------------------------------------
dldxx = spm_diff(funl,x,u,P,[1 1]);
for i = 1:l
    for j = 1:n
        D{i}(j,:) = dldxx{j}(i,:);
    end
end
for i = 1:l
    L2{i} = spm_cat(diag({0, D{i}}));
end
    

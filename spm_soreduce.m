function [M0,M1,M2,L1,L2] = spm_soreduce(M,P)
% reduction of a fully nonlinear MIMO system to second-order form
% FORMAT [M0,M1,M2,L1,L2] = spm_soreduce(M,P);
%
% M   - model specification structure
% Required fields:
%   M.f   - dx/dt    = f(x,u,P,M)        {function string or m-file}
%   M.g   - y(t)     = g(x,u,P,M)        {function string or m-file}
%   M.x   - (n x 1) = x(0) = expansion point: defaults to x = 0;
%   M.u   - (m x 1) = u    = expansion point: defaults to u = 0;
%
% P   - model parameters
%
% A second order approximation is returned where the states are
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
%       dq/dt = M0*q + 
%               u(1)*M1{1}*q + u(2)*M1{2}*q + ....
%               x(1)*M2{1}*q + x(2)*M2{2}*q + ....
%        y(i) = L(i,:)*q + ...
%
%--------------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_soreduce.m 5219 2013-01-29 17:07:07Z spm $


% set up
%==========================================================================

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

% add observer if not specified
%--------------------------------------------------------------------------
try
    fung = fcnchk(M.g,'x','u','P','M');
catch
    M.g  = inline('spm_vec(x)','x','u','P','M');
    M.l  = M.n;
    fung = fcnchk(M.g,'x','u','P','M');
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
[dfdxx,dfdx,f0] = spm_diff(funx,x,u,P,M,[1 1]);
[dfdxu,dfdx,f0] = spm_diff(funx,x,u,P,M,[1 2]);
dfdu  = spm_diff(funx,x,u,P,M,2);
m     = length(dfdxu);          % m inputs
n     = length(f0);             % n states



% Bilinear operators
%==========================================================================

% Bilinear operator - M0
%--------------------------------------------------------------------------
M0    = spm_cat({0                     []    ;
                (f0 - dfdx*spm_vec(x)) dfdx});

% Bilinear operator - M1 = dM0/du
%--------------------------------------------------------------------------
for i = 1:m
    M1{i} = spm_cat({0,                                []        ;
                    (dfdu(:,i) - dfdxu{i}*spm_vec(x)), dfdxu{i}});
end

% Bilinear operator - M2 = dM0/dx
%--------------------------------------------------------------------------
for i = 1:n
    M2{i} = spm_cat({0,                                []        ;
                    (dfdx(:,i) - dfdxx{i}*spm_vec(x)), dfdxx{i}});
end

if nargout < 4, return, end


% Output matrices - L1
%==========================================================================

% l(x(0),0)
%--------------------------------------------------------------------------
[dgdx,g0] = spm_diff(fung,x,u,P,M,1);
L1    = spm_cat({(g0 - dgdx*spm_vec(x)), dgdx});
l     = length(g0);

if nargout < 5, return, end

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
    

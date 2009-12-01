function [x] = spm_expm(J,x)
% approximate matrix exponential using a Taylor expansion
% FORMAT [y] = spm_expm(J,x)
% FORMAT [y] = spm_expm(J)
% y          = expm(J)*x:
% y          = expm(J);
%
% This routine covers and extends expm  functionality  by  using  a
% comoutationally  expedient  approximation  that can handle sparse
% matrices when dealing with the special case of expm(J)*x, where x
% is a vector, in an efficient fashion
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_expm.m 3605 2009-12-01 13:29:43Z karl $


% expm(J) use Pade approximation
%==========================================================================

% ensure norm is < 1/2 by scaling by power of 2
%--------------------------------------------------------------------------
J     = full(J);
I     = eye(size(J));
[f,e] = log2(norm(J,'inf'));
s     = max(0,e+1);
J     = J/2^s;
X     = J;
c     = 1/2;
E     = I + c*J;
D     = I - c*J;
q     = 16;
p     = 1;
for k = 2:q
    c   = c*(q - k + 1)/(k*(2*q - k + 1));
    X   = J*X;
    cX  = c*X;
    E   = E + cX;
    if p
        D = D + cX;
    else
        D = D - cX;
    end
    if norm(cX,1) < 1e-16;
        break
    end
    p = ~p;
end
E = D\E;  % E = inv(D)*E

% Undo scaling by repeated squaring E = E^(2^s)
%--------------------------------------------------------------------------
for k = 1:s
    E = E*E;
end

% Multiply by x if necessary
%--------------------------------------------------------------------------
if nargin > 1
    x = E*x;
else
    x = E;
end

return


% notes: Alterntive but slower evaluation of exp(J)*x
%==========================================================================

if nargin > 1

    % compute y = expm(J)*x = (1 + J + J*J/2! + J*J*J/3!  + ...)*x
    %----------------------------------------------------------------------
    J     = sparse(J);
    x     = sparse(x);
    x0    = x;
    fx    = J*x;
    j     = 1;
    nx    = norm(fx,1);
    TOL   = nx/256;

    while nx > TOL

        % accumulate high-order terms until convergence
        %------------------------------------------------------------------
        j  = j + 1;
        x  = x + fx;
        fx = J*fx/j;

        % revert to Pade approximation if numerical overflow
        %------------------------------------------------------------------
        if isinf(nx)
            fprintf('Reverting to Pade approximation')
            x  = spm_expm(J)*x0;
            return
        else
            nx = norm(fx,1);
        end
    end
    
else
    x  = spm_expm(J);
end

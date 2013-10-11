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
% $Id: spm_expm.m 5691 2013-10-11 16:53:00Z karl $

% % standard [eigen]solution
% %--------------------------------------------------------------------------
% [V,D] = eig(full(J));
% E     = V*diag(exp(diag(D)))/V;
%
% % Multiply by x if necessary
% %--------------------------------------------------------------------------
% if nargin > 1, x = E*x; else, x = E; end


% expm(J) use Pade approximation
%======================================================================

% ensure norm is < 1/2 by scaling by power of 2
%----------------------------------------------------------------------
I     = speye(size(J));
[f,e] = log2(norm(J,'inf'));
s     = max(0,e + 1);
J     = J/2^s;
X     = J;
c     = 1/2;
E     = I + c*J;
D     = I - c*J;
q     = 6;
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
    p = ~p;
end

% E = inv(D)*E
%--------------------------------------------------------------------------
E = D\E;  

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




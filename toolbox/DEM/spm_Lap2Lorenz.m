function [s,q,f] = spm_Lap2Lorenz(P,w,x)
% Laplace version of the Lorentz system
% FORMAT [s,q,f] = spm_Lap2Lorenz(P,[w,x])
% s - second order polynomial coefficients for(negative) potential
% q - second-order polynomial coefficients for flow
%
% this routine evaluates the Laplacian version of a Lorentz system with
% supplied parameters in terms of second order polynomial coefficients.
% This is an exact solution that conforms to the Helmholtz decomposition;
% however, with an improper steady-state density due to the absence of a
% leading diagonal Hessian. In the SPM code, the polynomial coefficients
% for the flow operator include coefficients for the leading diagonal. This
% means one can supplement any supplied dissipative or diagonal flow
% operator with state-dependent terms (e.g., state-dependent random
% fluctuations). The fixed values of these are specified in terms of the
% precision of random fluctuations (i.e., G = 1/(2*w)
%--------------------------------------------------------------------------
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% initialise polynomial coefficients
%--------------------------------------------------------------------------
s  = zeros(10,1,'like',P);
q  = zeros(60,1,'like',P);

% analytic form for potential (polynomial) coefficients
%--------------------------------------------------------------------------
q(44) = ((1 + P(1)) - sqrt((1 + P(1))^2 - 4*P(1)*P(2))) / (2*(P(2)/(1 + P(1)) - 1));
G(2)  = P(3)*(1 - q(44))/2;
G(1)  = G(2)*P(1)/P(3);

s(5)  = P(3)/G(2);
s(7)  = (1 + P(1))/q(44);

q(11) = G(1);
q(37) = 1/s(5);
q(55) = (q(44)*s(5) - 1)/s(7);
q(57) = P(2)/s(7);

G(3)  = (q(44) - q(57))/s(7);

% switch sign of flow coefficients
%--------------------------------------------------------------------------
q(1)  = -G(1);
q(31) = -G(2);
q(51) = -G(3);


return

% Notes for functional form of Laplace coefficients
%==========================================================================
sympref('FloatingPointOutput',1);
syms  x [1 3] 'real'
P     = [10; 8/3; 32];                  % parameters [sig, beta, rho]
[s,q] = spm_Lap2Lorenz(P);

% undo sign switch of flow coefficients
%--------------------------------------------------------------------------
G(1) = -q(1);
G(2) = -q(31);
G(3) = -q(51);

% evaluate flow and Helmholtz decomposition
%--------------------------------------------------------------------------
Q  = [ -G(1),               q(11),          0;
      -q(11), (q(37)*x(3) - G(2)), q(44)*x(2);
           0,         -q(44)*x(2), (q(57)*x(3) + q(55)*x(1)*x(2) - G(3))];
     
dS = [-s(5)*x(2);
      -s(5)*x(1);
          -s(7)];
     
L  = [0;
      0;
     q(44) - q(57)];
         
f = Q*dS - L

return






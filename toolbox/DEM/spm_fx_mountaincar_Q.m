function f = spm_fx_mountaincar_Q(x,v,P)
% state equations based on the Helmholtz decomposition
% FORMAT f = spm_fx_mountaincar_Q(x,v,P)
% x    - [x, x']
% v    - exogenous force
%
% P.a  - 0th order coefficients of Q
% P.b  - 1st order coefficients of Q
% P.c  - 2nd order coefficients of Q
%
% M    - model structure
%
% f    - flow dx/dt
%
% see:
% Gaussian Processes in Reinforcement Learning
% Carl Edward Rasmussen and Malte Kuss
% Max Planck Institute for Biological Cybernetics
% Spemannstraße 38, 72076 Tubingen, Germany
% {carl,malte.kuss}@tuebingen.mpg.de
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_mountaincar_Q.m 7679 2019-10-24 15:54:07Z spm $
 
% f = (D + Q(x))*grad(V)
%==========================================================================
 
% Q
%--------------------------------------------------------------------------
Q     = exp(-2) + P.Q(1) + P.Q(2)*x(1) + P.Q(3)*x(2) + P.Q(4)*x(1)*x(1) + P.Q(5)*x(1)*x(2) + P.Q(6)*x(2)*x(2);
Q     = [0 -Q ; 
         Q  0];
     
% grad(V)
%--------------------------------------------------------------------------
d     =  x - [1; 0];
p     =  P.V(1) + P.V(2)*x(1) + P.V(3)*x(2) + P.V(4)*x(1)*x(1) + P.V(5)*x(1)*x(2) + P.V(6)*x(2)*x(2);
Gp    = [P.V(2) + P.V(5)*x(2) + 2*P.V(4)*x(1);
         P.V(3) + P.V(5)*x(1) + 2*P.V(6)*x(2)];
J     = diag(exp(P.J));
V     = -(1/2)*(d'*d)*exp(p);
GV    = -exp(p)*d + V*Gp;
 
% flow (f)
%--------------------------------------------------------------------------
f     = J*GV + Q*GV;
 
% exogenous forces
%--------------------------------------------------------------------------
f(2)  = f(2) + v(1);

function [y] = spm_fx_fmri_dcm(x,u,P,M)
% state equation for a two state DCM of fMRI
% responses
% FORMAT [y] = spm_fx_dcm(x,u,P,M)
% x      - state vector
%   x(:,1) - excitatory neuronal activity             ue
%   x(:,2) - vascular signal                          s
%   x(:,3) - rCBF                                  ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoyxHb                               ln(q)
%   x(:,6) - inhibitory neuronal activity             ui
%
% y      - dx/dt
%__________________________________________________________________________
%
% References for state equations:
% 1. Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for 
%    fMRI: a two-state model.
%    Neuroimage. 2008 Jan 1;39(1):269-78. 
%
% 2. Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
%    Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
%    Neuroimage 42:649-662, 2008.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_fmri_dcm.m 3705 2010-02-01 20:51:28Z karl $
 
 
% hemodynamic parameters
%--------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) – auto-regulation                                d(ds/dt)/df)
%   H(3) - transit time                                   (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal

 
% neuronal states
%==========================================================================
 
% excitatory connections
%--------------------------------------------------------------------------
for i = 1:size(P.B,3)
    P.A = P.A + u(i)*P.B(:,:,i);
end
 
% and nonlinear (state) terms
%--------------------------------------------------------------------------
for i = 1:size(P.D,3)
    P.A = P.A + x(i,1)*P.D(:,:,i);
end
 
% extrinsic
%--------------------------------------------------------------------------
A     = exp(P.A)/8;             % enforce positivity
IE    = diag(diag(A));          % inhibitory to excitatory
EE    = A - IE;                 % excitatory to excitatory
EI    = 1;                      % excitatory to inhibitory
SE    = 1;                      % self-inhibition (excitatory)
SI    = 2;                      % self-inhibition (inhibitory)
 
% motion - excitatory and inhibitory: y = dx/dt
%--------------------------------------------------------------------------
y      = x;
y(:,1) = EE*x(:,1) - SE*x(:,1) - IE*x(:,6) + P.C*u;
y(:,6) = EI*x(:,1) - SI*x(:,6);
 
 
% hemodynamic states
%==========================================================================
 
% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(:,3:5)  = exp(x(:,3:5));
 
% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv        = x(:,4).^(1./P.H(:,4));
 
% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff        = (1 - (1 - P.H(:,5)).^(1./x(:,3)))./P.H(:,5);
 
% implement differential state equation y = dx/dt
%--------------------------------------------------------------------------
y(:,2)    = x(:,1)  - P.H(:,1).*x(:,2) - P.H(:,2).*(x(:,3) - 1);
y(:,3)    = x(:,2)./x(:,3);
y(:,4)    = (x(:,3) - fv)./(P.H(:,3).*x(:,4));
y(:,5)    = (ff.*x(:,3) - fv.*x(:,5)./x(:,4))./(P.H(:,3).*x(:,5));
y         = y(:);
 
 
return

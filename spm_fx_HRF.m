function [f] = spm_fx_HRF(x,u,P)
% state equation for the hemodynamic model
% FORMAT [f] = spm_fx_HRF(x,u,P)
% x      - state vector
%   x(1) - vascular signal           (s)
%   x(2) - rCBF                      (f)
%   x(3) - venous volume             (v)
%   x(4) - deoyxHb                   (q)
% u      - input (neuronal activity) (u)
% P      - free parameter vector
%   P(1) - signal decay      - d(ds/dt)/ds)
%   P(2) - autoregulation    - d(ds/dt)/df)
%   P(3) - transit time                (t0)
%   P(4) - exponent for Fout(v)     (alpha)
%   P(5) - resting oxygen extraction   (E0)
%
%   P(5 + 1:m)   - input efficacies  - d(ds/dt)/du)
%
% y      - dx/dt
%__________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_fx_HRF.m 455 2006-02-22 18:40:33Z karl $


% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv     = x(3)^(1/P(4));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff     = (1 - (1 - P(5))^(1/x(2)))/P(5);

% implement differential state equation
%--------------------------------------------------------------------------
f(1)   = P(6:end)'*u - P(1)*x(1) - P(2)*(x(2) - 1);
f(2)   = x(1);
f(3)   = (x(2) - fv)/P(3);
f(4)   = (ff*x(2) - fv*x(4)/x(3))/P(3);
f      = f(:);

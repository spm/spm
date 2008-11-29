function [y] = spm_fx_dcm(x,u,P,M)
% state equation for a dynamic [bilinear/nonlinear/Balloon] model of fMRI 
% responses
% FORMAT [y] = spm_fx_dcm(x,u,P,M)
% x      - state vector
%   x(:,1) - neuronal acivity                         u
%   x(:,2) - vascular signal                          s
%   x(:,3) - rCBF                                  ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoyxHb                               ln(q)
% y      - dx/dt
%___________________________________________________________________________
%
% References for hemodynamic & neuronal state equations:
% 1. Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864, 1998.
% 2. Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in 
% fMRI: the Balloon model, Volterra kernels, and other hemodynamics. 
% Neuroimage 12:466-477, 2000.
% Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE, Breakspear 
% M, Friston KJ. Nonlinear dynamic causal models for fMRI. Neuroimage 42:
% 649-662, 2008. 
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Klaas Enno Stephan
% $Id: spm_fx_dcm.m 2504 2008-11-29 15:53:11Z klaas $


% hemodynamic parameters
%---------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) - autoregulation                                 d(ds/dt)/df)
%   H(3) - transit time                                   (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal   

% get dimensions
%---------------------------------------------------------------------------
if size(u,2) > size(u,1)                % make sure u is a column vector
    u = u';
end
if size(x,2) > size(x,1)                % make sure x is a column vector
    x = x';
end
m         = size(u,1);                  % number of inputs
n         = size(x,1)/5;                % number of regions

% reshape parameters
%---------------------------------------------------------------------------
if ~M.nlDCM
    % bilinear DCM
    [A B C H]   = spm_dcm_reshape(P,m,n);
else
    % nonlinear DCM
    [A B C H D] = spm_dcm_reshape(P,m,n);
end

% effective intrinsic connectivity
%---------------------------------------------------------------------------
for i = 1:m
      A = A + u(i)*B(:,:,i);
end

% configure state variables
%---------------------------------------------------------------------------
x         = reshape(x,n,5);

% modulatory (2nd order) terms
%---------------------------------------------------------------------------
if M.nlDCM
    % nonlinear DCM
    xD = zeros(size(A,1),size(A,1));
    for j = 1:n,
        xD = xD + x(j,1)*D(:,:,j);
    end
end

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(:,3:5)  = exp(x(:,3:5));  

% Fout = f(v) - outflow
%---------------------------------------------------------------------------
fv        = x(:,4).^(1./H(:,4));

% e = f(f) - oxygen extraction
%---------------------------------------------------------------------------
ff        = (1 - (1 - H(:,5)).^(1./x(:,3)))./H(:,5);

% implement differential state equation y = dx/dt
%---------------------------------------------------------------------------
y         = zeros(n,5);
if ~M.nlDCM
    % bilinear DCM
    y(:,1)    = A*x(:,1) + C*u;
else
    % nonlinear DCM
    y(:,1)    = A*x(:,1) + C*u + xD*x(:,1);
end
y(:,2)    = x(:,1)  - H(:,1).*x(:,2) - H(:,2).*(x(:,3) - 1);
y(:,3)    = x(:,2)./x(:,3);
y(:,4)    = (x(:,3) - fv)./(H(:,3).*x(:,4));
y(:,5)    = (ff.*x(:,3) - fv.*x(:,5)./x(:,4))./(H(:,3).*x(:,5));
y         = y(:);

return

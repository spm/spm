function [y] = spm_fx_dcm(x,u,P)
% state equation for a dynamic [Bilinear/Balloon] model of fMRI responses
% FORMAT [y] = spm_fx_dcm(x,u,P)
% x      - statw vector
%   x(:,1) - neuronal acivity
%   x(:,2) - vascular signal           (s)
%   x(:,3) - rCBF                      (f)
%   x(:,4) - venous volume             (v)
%   x(:,5) - deoyxHb                   (q)
% y      - dx/dt
%___________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%___________________________________________________________________________
% %W% Karl Friston %E%


% hemodynamic parameters
%---------------------------------------------------------------------------
%   H(1) - signal decay      - d(ds/dt)/ds)
%   H(2) - autoregulation    - d(ds/dt)/df)
%   H(3) - transit time                (t0)
%   H(4) - exponent for Fout(v)     (alpha)
%   H(5) - resting oxygen extraction   (E0)

H         = [0.65 0.41 0.98 0.32 0.34];

% get dimensions
%---------------------------------------------------------------------------
m         = size(u,1);
r         = size(x,1)/5;

% connectivity
%---------------------------------------------------------------------------
[A B C D] = spm_bl_reshape(full(P),m,r,0);
for     i = 1:m
	A = A + u(i)*B(:,:,i);
end

% configure state variables
%---------------------------------------------------------------------------
x         = reshape(x,r,5);
x(:,3:5)  = x(:,3:5) + 1;

%    Fout = f(v) - outflow
%---------------------------------------------------------------------------
fv        = x(:,4).^(1/H(4));

% e = f(f) - oxygen extraction
%---------------------------------------------------------------------------
ff        = (1 - (1 - H(5)).^(1./x(:,3)))/H(5);

% implement differential state equation y = dx/dt
%---------------------------------------------------------------------------
y         = zeros(r,5);
y(:,1)    = A*x(:,1) + C*u + D;
y(:,2)    = x(:,1) - H(1)*x(:,2) - H(2)*(x(:,3) - 1);
y(:,3)    = x(:,2);
y(:,4)    = (x(:,3) - fv)/H(3);
y(:,5)    = (ff.*x(:,3) - fv.*x(:,5)./x(:,4))/H(3);
y         = y(:);

function f = spm_fx_mountaincar(x,v,varargin)
% state equations for mountain car problem
% FORMAT f = spm_fx_mountaincar(x,v,P)
% FORMAT f = spm_fx_mountaincar(x,v,a,P)
% FORMAT f = spm_fx_mountaincar(x,v,P,M)
% x   - [x, x']
% v   - exogenous force
% a   - action
%
% P.a - 0th order coefficients of force
% P.b - 1st order coefficients of force
% P.c ñ 2nd order coefficients of force
% P.d - action coefficient
%
% see:
% Gaussian Processes in Reinforcement Learning
% Carl Edward Rasmussen and Malte Kuss
% Max Planck Institute for Biological Cybernetics
% Spemannstraﬂe 38, 72076 T®ubingen, Germany
% {carl,malte.kuss}@tuebingen.mpg.de
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_mountaincar.m 3122 2009-05-13 16:05:59Z karl $
 
 
% determine controlled forces (a)
%==========================================================================
 
% spm_fx_mountaincar(x,v,P) - recognition model (no action)
%--------------------------------------------------------------------------
if length(varargin) == 1
    
    a    = 0;
    P    = varargin{1};
    
else
    % spm_fx_mountaincar(x,v,P,M) - generative model (no action)
    %----------------------------------------------------------------------
    try
        varargin{2}.f;
        varargin{2}.g;
        a    = 0;
        P    = varargin{1};
 
    % spm_fx_mountaincar(x,v,a,P) - generative model (action)
    %----------------------------------------------------------------------
    catch
        a    = varargin{1};
        P    = varargin{2};
    end
end
 
 
% acceleration = force:
%--------------------------------------------------------------------------
a     = tanh((P.d*a + P.a + P.b*x + P.c*kron(x,x))/2);
 
% f(x)
%--------------------------------------------------------------------------
dt    = 1/4;
if x(1) < 0                                  % gravity
    dHdx = 2*x(1) + 1;
else
    xx   = x(1)^2;
    dHdx = (1 + 5*xx)^(-1/2) - 5*xx/(1 + 5*xx)^(3/2) + (x(1)/2)^4;
end
f     = [x(2); a + v - dHdx - x(2)/8]*dt;
 
return
 
% NOTES: Plots for figure
%--------------------------------------------------------------------------
dx    = 1/64;
x     = linspace(-2,2,1/dx);
xx    = x.^2;
dHdx  = (x < 0).*(2*x + 1);
dHdx  = (x > 0).*(1./(1 + 5*xx).^(1/2) - 5*xx./(1 + 5*xx).^(3/2) + (x/2).^4) + dHdx;
H     = cumsum(dHdx)*dx;
H     = H - min(H);
 
subplot(2,2,1)
plot(x,H,x.^0,H,':')
xlabel('position','FontSize',12)
ylabel('height','FontSize',12)
title('mountain car problem','FontSize',16)
axis square
 
subplot(2,2,2)
plot(x,-dHdx,x,x*0,'-.',x,-x.^0,'--')
xlabel('position','FontSize',12)
ylabel('force','FontSize',12)
title('forces','FontSize',16)
axis square

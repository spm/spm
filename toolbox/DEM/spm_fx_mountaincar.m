function f = spm_fx_mountaincar(x,v,a,P)
% state equations for mountain car problem
% FORMAT f = spm_fx_mountaincar(x,v,P)
% FORMAT f = spm_fx_mountaincar(x,v,a,P)
% FORMAT f = spm_fx_mountaincar(x,v,P,M)
% x    - [x, x']
% v    - exogenous force
% a    - action
%
% P.a  - 0th order coefficients of force
% P.b  - 1st order coefficients of force
% P.c  - 2nd order coefficients of force
% P.d  - action coefficient
%
% M    - model structure
%
% f    - flow dx/dt
%
% see:
% Gaussian Processes in Reinforcement Learning
% Carl Edward Rasmussen and Malte Kuss
% Max Planck Institute for Biological Cybernetics
% Spemannstraße 38, 72076 T¨ubingen, Germany
% {carl,malte.kuss}@tuebingen.mpg.de
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_mountaincar.m 7679 2019-10-24 15:54:07Z spm $
 
 
% determine controlled forces (a)
%==========================================================================
 
% spm_fx_mountaincar(x,v,P) - recognition model (no action)
%--------------------------------------------------------------------------
global eta
if isempty(eta), eta = 8; end
if nargin == 3
    
    P    = a;
    a    = 0;
  
elseif all(isfield(P,{'f','g'}))
    
    % spm_fx_mountaincar(x,v,P,M) - generative model (no action)
    %----------------------------------------------------------------------
    P.f;
    P.g;
    M    = P;
    P    = a;
    a    = 0;

end


% default parameters
%--------------------------------------------------------------------------
if isempty(P)
    P.a  = 0;
    P.b  = [0 0];
    P.c  = [0 0 0 0];
    P.d  = 1;
end
 
% acceleration = force:
%--------------------------------------------------------------------------
x     = x(1:2);
a     = tanh((P.d*a + P.a + P.b*x + P.c*kron(x,x))/2);
 
% f(x)
%--------------------------------------------------------------------------
dt    = 1/4;
if x(1) < 0                                  % gravity
    dHdx = 2*x(1) + 1;
else
    xx   = x(1)^2;
    dHdx = (5*xx + 1).^(-3/2) + (x(1)/2).^4;
end
f     = [x(2); a + v - dHdx - x(2)/eta]*dt;
 

return

 
% NOTES: Plots for figure
%--------------------------------------------------------------------------
dx    = 1/64;
x     = linspace(-2,2,1/dx);
xx    = x.^2;
dHdx  = (x < 0).*(2*x + 1);
dHdx  = (x > 0).*((5*xx + 1).^(-3/2) + (x/2).^4) + dHdx;
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

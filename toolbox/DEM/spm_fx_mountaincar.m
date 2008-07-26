function f = spm_fx_mountaincar(x,v,varargin)
% state equations for mountain car problem
% FORMAT f = spm_fx_mountaincar(x,v,P)
% FORMAT f = spm_fx_mountaincar(x,v,a,P)
% FORMAT f = spm_fx_mountaincar(x,v,P,M)
% x   - [x, x']
% v   - ezogenous force
% a   - action
%
% P.a - 0th order coefficents of force
% P.b - 1st order coefficents of force
% P.c - 2nd order coefficents of force
% P.d - action coefficent
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
% $Id: spm_fx_mountaincar.m 1961 2008-07-26 09:38:46Z karl $


% get dimensions and configure state variables
%--------------------------------------------------------------------------
if length(varargin) == 1
    a = 0;
    P = varargin{1};
else
    % spm_fx_mountaincar(x,v,P,M)
    %----------------------------------------------------------------------
    try
        varargin{2}.f;
        varargin{2}.g;
        a = 0;
        P = varargin{1};
        
    % spm_fx_mountaincar(x,v,a,P)
    %----------------------------------------------------------------------
    catch
        a = varargin{1};
        P = varargin{2};
    end
end

% acceleration = force: |a| < 4
%--------------------------------------------------------------------------
a     = 2*spm_phi(P.a + P.b*x + P.c*spm_kron(x,x) + P.d*a) - 1;

% f(x)
%--------------------------------------------------------------------------
dt    = 1/4;
if x(1) < 0                                  % gravity
    dHdx = 2*x(1) + 1;
else
    xx   = x(1)^2;
    dHdx = 1/(1 + 5*xx)^(1/2) - 5*xx/(1 + 5*xx)^(3/2);
end

f     = [x(2);
         a + v - dHdx - x(2)/4 - x(1)/4]*dt;
    

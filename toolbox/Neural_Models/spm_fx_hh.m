function [y] = spm_fx_hh(x,u,P)
% state equation for a single Hodgkin-Huxley like unit
% FORMAT [y] = spm_fx_hh(x,u,P)
%
% states
%--------------------------------------------------------------------------
% x(1) = proportion of open channels        % AMPA
% x(2) = proportion of open channels        % GABA
% x(3) = proportion of open channels        % K - slow
% x(4) = proportion of open channels        % NMDA
% x(5) = V                  % transmembrane potential mV
% x(6) = t                  % time since last spike
%
% u    = input - opening rate of AMPA channels
%
% P(1) = opening rate of AMPA channels
% P(1) = opening rate of GABA channels
% P(1) = opening rate of NMDA channels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_hh.m 5615 2013-08-15 14:37:24Z spm $


% fixed parameters
%--------------------------------------------------------------------------
C     =   0.375;                % Capacitance {nF}
Vl    =     -73;                % Resting potential {mV}
gl    =      25;                % passive conductance {nS}
Vd    =     -53;                % depolarization potential {mV}

% channel: AMPA GABA   K-s  NMDA
%--------------------------------------------------------------------------
V     = [   00  -70   -90  00];         % Equilibrium potential {mV}
g     = [   24   64   128   8];         % active  conductance {nS}
t     = [  2.4    7    80  100]*1e-3;       % time constant {ms}


% dV/dt {mV/s)
%--------------------------------------------------------------------------
v     = x(5);
T     = x(6);
dVdt  = (1/C) * (gl*(Vl   - v) + ...                % leak
     g(1)*x(1).*(V(1) - v) + ...                % AMPA
     g(2)*x(2).*(V(2) - v) + ...                % GABA
     g(3)*x(3).*(V(3) - v) + ...                % K-slow
     g(4)*x(4).*(V(4) - v)/(1 + exp(-(v + 10)/14)) );   % NMDA

s     = exp(-T^2/(2*(1e-3)^2));
dVdt  = dVdt + 1e4*(-90 - v)*s;

% pst
%--------------------------------------------------------------------------
dtdt  = 1 - T*1e4*(v > Vd);


% dp/dt, u-induced opening - spontaneous closing
%--------------------------------------------------------------------------
dpedt = ((1 - x(1))*(P(1) + u) - x(1))/t(1);
dpidt = ((1 - x(2))*(P(2)    ) - x(2))/t(2);
dpkdt = ((1 - x(3))*(4*s     ) - x(3))/t(3);
dpndt = ((1 - x(4))*(P(3)    ) - x(4))/t(4);


% dx/dt
%--------------------------------------------------------------------------
y     = [dpedt; dpidt; dpkdt; dpndt; dVdt; dtdt];

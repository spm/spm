function [y] = spm_gx_hh(x,u,P)
% output for a single Hodgkin-Huxley like unit
% FORMAT [y] = spm_gx_hh(x,u,P)
%
% outputs
%--------------------------------------------------------------------------
% y(1) = V transmembrane potential mV (c.f. LFP)
% y(2) = spike rate (Hz) = 1/PST
% y(3) = dendritic energy = g(1)*x(1).*(V(1) - v).^2 + ... (mV.^2mS)
%--------------------------------------------------------------------------
%
% states
%--------------------------------------------------------------------------
% x(1) = proportion of open channels        % AMPA
% x(2) = proportion of open channels        % GABA
% x(3) = proportion of open channels        % K - slow
% x(4) = proportion of open channels        % NMDA
% x(5) = V                  % transmembrane potential mV
% x(6) = t                  % time since last spike (peri-spike time)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_hh.m 4936 2012-09-18 19:47:55Z karl $
 
 
% fixed parameters
%--------------------------------------------------------------------------
Vl    =     -73;                % Resting potential {mV}
gl    =      25;                % passive conductance {nS}
 
% channel: AMPA GABA   K-s  NMDA
%--------------------------------------------------------------------------
V     = [   00  -70   -90  00];         % Equilibrium potential {mV}
g     = [   24   64   128   8];         % active  conductance {nS}
 
% energy = g(i)*x(8)*(V(i) - v)
%--------------------------------------------------------------------------
v     = x(5);
T     = x(6);
s     = exp(-T^2/(2*(1e-3)^2));
E     =  gl*(Vl   - v).^2 + ...                 % leak
     g(1)*x(1).*(V(1) - v).^2 + ...             % AMPA
     g(2)*x(2).*(V(2) - v).^2 + ...             % GABA
     g(3)*x(3).*(V(3) - v).^2 + ...             % K-slow
     g(4)*x(4).*(V(4) - v).^2/(1 + exp(-(v + 10)/14));  % NMDA
 
% output (y)
%--------------------------------------------------------------------------
y      = [v; 1/(1/128 + T); E/1000];



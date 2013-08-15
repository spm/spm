function s = mg_switch(V)
% switching output s: determined by voltage (V) depedent magnesium blockade
% parameters as per Durstewitz, Seamans & Sejnowski 2000

% $Id: mg_switch.m 5615 2013-08-15 14:37:24Z spm $

s = 1.50265./(1 + 0.33*exp(-0.06.*V));

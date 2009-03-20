function [G] = spm_lx_phase (P,M)
% Observation function for phase-coupled oscillators
% FORMAT [G] = spm_lx_phase (P,M)
%
% G     Observations y = Gx

Nr=length(P.L);
G=eye(Nr);
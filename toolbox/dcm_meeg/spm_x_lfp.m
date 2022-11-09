function [x] = spm_x_lfp(P)
% Initial state of a neural mass model of erps
% FORMAT [x] = spm_x_lfp(P)
% P - parameters
%
% x        - x(0)
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% array of states
%--------------------------------------------------------------------------
n  = length(P.A{1});                          % number of sources
m  = 13;                                      % number of states
x  = sparse(n,m);

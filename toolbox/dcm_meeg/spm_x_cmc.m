function [x] = spm_x_cmc(P)
% returns the initial state of a canonical microcircuit model
% FORMAT [x] = spm_x_cmc(P)
% P - parameters
%
% x        - x(0)
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_x_cmc.m 4232 2011-03-07 21:01:16Z karl $
 
% array of states
%--------------------------------------------------------------------------
n  = length(P.A{1});                          % number of sources
m  = 8;                                       % number of states
x  = sparse(n,m);


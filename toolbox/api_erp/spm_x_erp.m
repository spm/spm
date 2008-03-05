function [x] = spm_x_erp(P)
% returns the initial state of a neural mass model of erps
% FORMAT [x] = spm_x_erp(P)
% M - model structure
% P - parameters
%
% x        - x(0)
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_x_erp.m 1189 2008-03-05 17:19:26Z karl $

% array of states
%--------------------------------------------------------------------------
n  = length(P.A{1});                          % number of sources
m  = 9;                                       % number of states
x  = sparse(n,m);


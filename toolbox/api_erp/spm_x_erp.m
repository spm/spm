function [x] = spm_x_erp(P,M,U)
% returns the initial state of a neural mass model of erps
% FORMAT [x] = spm_x_erp(P,M,U)
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
% $Id: spm_x_erp.m 1174 2008-02-27 20:22:30Z karl $

% P.x0 = spiny stellate depolarisation
%--------------------------------------------------------------------------
n     = length(P.A{1});                       % number of sources
m     = 9;                                    % number of states
try
    x = sparse([1:n],1,P.x0,n,m);
catch
    x = sparse(n,m);
end

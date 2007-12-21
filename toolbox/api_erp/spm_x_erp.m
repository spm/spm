function [x] = spm_x_erp(P,M,U)
% reutrn the intial state of a neural mass model of erps
% FORMAT [x] = spm_x_erp(P,M,U)
% M - model structure
% P - parameters
%
% x        - x(0)
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

% spiny stellate depolarisation
%--------------------------------------------------------------------------
n     = length(P.A{1});
try
    x = sparse([1:n] + 1,1,P.x0,n*9 + 1,1);
catch
    x = sparse(n*9 + 1,1);
end

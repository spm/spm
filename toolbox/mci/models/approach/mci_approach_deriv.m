function [dLdp,iCpY,L] = mci_approach_deriv (P,M,U,Y)
% Gradient of log-likelihood for approach model
% FORMAT [dLdp,iCpY,L] = mci_approach_deriv (P,M,U,Y)
%
% dLdp      gradient of log joint
% iCpY      curvature (Fisher Information)
% L         log joint
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_approach_deriv.m 6548 2015-09-11 12:39:47Z will $

G = mci_approach_gen (P,M,U);
if isstruct(Y)
    e = Y.y-G;
else
    e = Y-G;
end

V=exp(P(1));
tau=exp(P(2));
t=U.X;

y=-60+V*(1-exp(-t/tau));

dydp = [V*(1-exp(-t/tau)),-V*exp(-t/tau).*(t/tau)];
dLdp = dydp'*M.iCe*e;
iCpY = dydp'*M.iCe*dydp;

if nargout > 2
    L = mci_approach_like (P,M,U,Y);
end



function [dLdp,iCpY,L] = mci_linear_deriv (P,M,U,Y)
% Gradient of likelihood for linear regression
% FORMAT [dLdp,iCpY,L] = mci_linear_deriv (P,M,U,Y)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linear_deriv.m 6275 2014-12-01 08:41:18Z will $

G = mci_linear_gen (P,M,U);
if isstruct(Y)
    e = Y.y-G;
else
    e = Y-G;
end
X = U.X;

dLdp=X'*M.iCe*e;
dLdp=dLdp';

iCpY=X'*M.iCe*X;

if nargout > 2
    L=mci_linear_like (P,M,U,Y);
end
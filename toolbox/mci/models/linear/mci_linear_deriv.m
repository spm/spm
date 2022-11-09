function [dLdp,iCpY,L] = mci_linear_deriv (P,M,U,Y)
% Gradient of likelihood for linear regression
% FORMAT [dLdp,iCpY,L] = mci_linear_deriv (P,M,U,Y)
%
% P         parameters
% M         model
% U         inputs
% Y         data
%
% dLdp      gradient of log joint
% iCpY      curvature (Fisher Information)
% L         log joint
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

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
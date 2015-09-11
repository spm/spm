function [dLdp,iCpY,L] = mci_linsqr_deriv (P,M,U,Y)
% Gradient of likelihood for linear regression
% FORMAT [dLdp,iCpY,L] = mci_linsqr_deriv (P,M,U,Y)
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
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linsqr_deriv.m 6548 2015-09-11 12:39:47Z will $

G = mci_linsqr_gen (P,M,U);
if isstruct(Y)
    e = Y.y-G;
else
    e = Y-G;
end
X = U.X;

N=size(X,1);
dydp=2*(ones(N,1)*P').*X;

%dydp = spm_diff(M.IS,P,M,U,1);

dLdp=dydp'*M.iCe*e;
dLdp=dLdp';

iCpY=dydp'*M.iCe*dydp;

if nargout > 2
    L=mci_linsqr_like (P,M,U,Y);
end
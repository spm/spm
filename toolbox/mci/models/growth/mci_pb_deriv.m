function [dLdp,iCpY,L] = mci_pb_deriv (P,M,U,Y)
% Gradient of log-likelihood for Preece-Baines model
% FORMAT [dLdp,iCpY,L] = mci_pb_deriv (P,M,U,Y)
%
% dLdp      gradient of log joint
% iCpY      curvature (Fisher Information)
% L         log joint
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_pb_deriv.m 6548 2015-09-11 12:39:47Z will $

dydp = spm_diff(M.IS,P,M,U,1);

G = mci_pb_gen (P,M,U);
if isstruct(Y)
    e = Y.y-G;
else
    e = Y-G;
end

N=length(e);
iCe=M.iCe*eye(N);

dLdp = dydp'*iCe*e;
%dLdp = spm_diff(M.L,P,M,U,Y,1);

iCpY = dydp'*iCe*dydp;

if nargout > 2
    L = mci_pb_like (P,M,U,Y);
end



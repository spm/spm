function [dLdp,iCpY,L] = mci_discount_deriv (P,M,U,Y)
% Gradient of likelihood for discounting model
% FORMAT [dLdp,iCpY,L] = mci_discount_deriv (P,M,U,Y)
%
% P         parameters
% M         model structure
% U         contains rewards and times
% Y         data
%
% dLdp      gradient of log joint
% iCpY      curvature (Fisher Information)
% L         log joint
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_discount_deriv.m 6548 2015-09-11 12:39:47Z will $

dLdp = spm_diff(M.L,P,M,U,Y,1);
dLdp = full(dLdp(:));

X = spm_diff('mci_discount_act',P,M,U,1);
g = mci_discount_gen (P,M,U);
Lambda=diag(g.*(1-g));
iCpY=X'*Lambda*X;

if nargout > 2
    L=mci_discount_like (P,M,U,Y);
end
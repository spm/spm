function [dLdp,iCpY,L] = spm_mci_diff (P,M,U,Y)
% Compute gradient and curvature of log likelihood using finite differences
% FORMAT [dLdp,iCpY,L] = spm_mci_diff (P,M,U,Y)
%
% dLdp      gradient of log likelihood
% iCpY      curvature (observed Fisher Information)
% L         log likelihood
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_diff.m 6548 2015-09-11 12:39:47Z will $

dLdp = spm_diff(M.L,P,M,U,Y,1);
dLdp = full(dLdp(:))';

if nargout > 1
    % Compute Hessian
    H=spm_diff(M.L,P,M,U,Y,[1 1]);
    Np=length(H);
    for p=1:Np,
        iCpY(:,p)=-full(H{p}(:));
    end
end

if nargout > 2
    L = feval(M.L,P,M,U,Y);
end

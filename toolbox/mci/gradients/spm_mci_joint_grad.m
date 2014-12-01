function [j,iCpY,st,L] = spm_mci_joint_grad (Pr,M,U,Y)
% Gradient of Log Joint Probability
% FORMAT [j,iCpY,st,L] = spm_mci_joint_grad (Pr,M,U,Y)
%
% Pr        parameters (vectorised and in M.V subspace)
% M         model structure
% U         inputs
% Y         data
%
% j         gradient of log joint, dL/dP 
% iCpY      Inverse Cov of params due to noise on outputs
% st        Status flag (0 for OK, -1 for problem)
% L         log joint, L = log p(Y,P)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_joint_grad.m 6275 2014-12-01 08:41:18Z will $

st=0;
Pr=Pr(:);

dLdp_prior = spm_mci_gprior_deriv (Pr,M);

% Parameters in original space
P = M.V*Pr+M.vpE;

% g     gradient of log likelihood, d[log p(Y|P)]/dP
if nargout > 1
    if isfield(M,'dL')
        [g,iCpY,L2] = feval(M.dL,P,M,U,Y);
        % Project into eigenparam space
        g=M.V'*g(:);
        g=g';
        iCpY=M.V'*iCpY*M.V;
    else
        [g,iCpY,st,L2] = spm_mci_glike_deriv (P,M,U,Y);
    end
else
    if isfield(M,'dL')
        g = feval(M.dL,P,M,U,Y);
        % Project into eigenparam space
        g=M.V'*g(:);
        g=g';
    else
        g = spm_mci_glike_deriv (P,M,U,Y);
    end
end

j = dLdp_prior + g;

if nargout > 3
    L1 =- e'*M.ipC*e/2 + M.log_prior_t2;
    L=L1+L2;
end

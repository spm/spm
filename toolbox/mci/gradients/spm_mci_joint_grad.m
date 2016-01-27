function [j,iCpY,st,L,L2] = spm_mci_joint_grad (Pr,M,U,Y)
% Gradient of Log Joint Probability
% FORMAT [j,iCpY,st,L,L2] = spm_mci_joint_grad (Pr,M,U,Y)
%
% Pr        parameters (vectorised and in M.V subspace)
% M         model structure. If field .beta is specified this
%           sets the inverse temperature to beta (default=1)
% U         inputs
% Y         data
%
% j         gradient of log joint, dL/dP 
% iCpY      Curvature (Fisher Information)
% st        Status flag (0 for OK, -1 for problem)
% L         log joint, L = log p(Y,P)
% L2        log likelihood, L2 = log p(Y|P)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_joint_grad.m 6697 2016-01-27 14:57:28Z spm $

try, beta=M.beta; catch, beta=1; end
    
st=0;
Pr=Pr(:);

dLdp_prior = spm_mci_gprior_deriv (Pr,M);

% Parameters in original space
P = M.V*Pr+M.vpE;

% g     gradient of log likelihood, d[log p(Y|P)]/dP
if nargout > 1
    if isfield(M,'dL')
        % User specified routines
        [g,iCpY,L2] = feval(M.dL,P,M,U,Y);
        % Project into eigenparam space
        g=M.V'*g(:);
        g=g';
        iCpY=M.V'*iCpY*M.V;
    else
        if isfield(M,'f')
            % For dynamical models use forward sensitivity approach
            [g,iCpY,st,L2] = spm_mci_glike_deriv (P,M,U,Y);
        else
            % For other models use finite differences
            [g,iCpY,L2] = spm_mci_diff(P,M,U,Y);
        end
    end
    iCpY=beta*iCpY;
else
    if isfield(M,'dL')
        % User specified routines
        g = feval(M.dL,P,M,U,Y);
        % Project into eigenparam space
        g=M.V'*g(:);
        g=g';
    else
        if isfield(M,'f')
            % For dynamical models use forward sensitivity approach
            g = spm_mci_glike_deriv (P,M,U,Y);
        else
            % For other models use finite differences
            g = spm_mci_diff(P,M,U,Y);
        end
    end
end

j = dLdp_prior + beta*g;

if nargout > 3
    e = Pr;
    L1 = - e'*M.ipC*e/2 + M.log_prior_t2;
    L = L1+beta*L2;
end

function [dLdp,g,x] = spm_mci_adjoint (Pr,M,U,Y)
% Gradient of log joint from adjoint method 
% FORMAT [dLdp,g,x] = spm_mci_adjoint (Pr,M,U,Y)
%
% Pr        Parameters (vectorised and in M.V subspace)
% M         Model structure
% U         Inputs  [Nin x N]
% Y         Data
%     
% dLdp      Gradient    [Np x 1]
% g         Outputs     [N x Nout]
% x         States      [N x Nstates]
%
% If M.adjlike=1 this function returns gradient of log likelihood
%
% This function uses integrators from MATLAB's ODE Suite
%
% B. Sengupta, K. Friston and W. Penny (2014) Efficient Gradient
% Computation for Dynamical Models. Neuroimage,98, 521-527. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_adjoint.m 6697 2016-01-27 14:57:28Z spm $

try, adjlike=M.adjlike; catch, adjlike=0; end

% Parameters in original space
P = M.V*Pr+M.vpE;

if isempty(U)
    U=zeros(1,M.N);
end

[g,x] = spm_mci_fwd(P,M,U);

% Gradient at computed times
% When computing output sensitivities, assume dydx=L, ie not a
% function of x. Generalise later
[tmp,L]=feval(M.g,M.x0,U(:,1),P,M);
e=Y-g;
djdx=e*M.iCe*L;

% integrate adjoint equation
lambda = spm_mci_adjoint_int(U,P,M,x,djdx);

% If observation function becomes dependent on parameters
% djdp will need updating
djdp = zeros(1,M.Np);

dLdp=zeros(1,M.Np);
for n=1:M.N,
    % Evaluate parameter Jacobian
    if isfield(M,'dfdp')
        Fp = feval(M.dfdp,x(n,:)',U(:,n),P,M);
    else
        Fp = spm_diff(M.f,x(n,:)',U(:,n),P,M,3);
    end
    dLdp=dLdp+djdp-lambda(n,:)*Fp;
end

if ~adjlike
    dlogpriordp = spm_mci_gprior_deriv (Pr,M);
    dLdp=dLdp+dlogpriordp;
end


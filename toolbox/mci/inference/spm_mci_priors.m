function [M] = spm_mci_priors (M)
% Quantities for computing log prior in subspace
% FORMAT [M] = spm_mci_priors (M)
%
% M.V               projection matrix
% M.ipC             Inverse prior cov in reduced space
% M.log_prior_t2    second term of log prior 
% M.Np              dimension of reduced space
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

if isstruct(M.pC)
    pC=full(diag(spm_vec(M.pC)));
else
    pC = M.pC;
end
V  = spm_svd(pC,exp(-32));
Np = size(V,2);
pC = V'*pC*V;
ipC = inv(pC);
log_prior_t2 = spm_logdet(ipC)/2-0.5*Np*log(2*pi);

M.ipC=ipC;
M.V=V;
M.log_prior_t2=log_prior_t2;
M.Np=Np;

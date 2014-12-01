function [M] = spm_mci_switch_prep (M)
% Prepare quantities for computing log prior in SVD-reduced space
% FORMAT [M] = spm_mci_switch_prep (M)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_switch_prep.m 6275 2014-12-01 08:41:18Z will $

% Big model

pC = M{1}.pC;
Np = size(pC,1);
V  = spm_svd(pC,exp(-32));
pC = V'*pC*V;
ipC = inv(pC);
log_prior_t2 = spm_logdet(ipC)/2-0.5*Np*log(2*pi);

M{1}.ipC=ipC;
M{1}.V=V;
M{1}.log_prior_t2=log_prior_t2;

% Small model

pC = M{2}.pC;
pC = V'*pC*V;
ipC = inv(pC);

ind=M{2}.subset;
p2=length(ind);
log_prior_t2 = spm_logdet(ipC(ind,ind))/2-0.5*p2*log(2*pi);
%log_prior_t2 = spm_logdet(ipC)/2-0.5*p2*log(2*pi);

M{2}.ipC=ipC;
M{2}.log_prior_t2=log_prior_t2;
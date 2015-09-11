function [M] = spm_mci_switch_prep (M)
% Prepare quantities for computing log prior in SVD-reduced space
% FORMAT [M] = spm_mci_switch_prep (M)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_switch_prep.m 6548 2015-09-11 12:39:47Z will $

for i=1:2,
    if isstruct(M{i}.pE)
        M{i}.vpE=spm_vec(M{i}.pE);
    else
        M{i}.vpE=M{i}.pE;
    end
    if isfield(M{i},'Ce')
        M{i}.logdet_Ce=spm_logdet(M{i}.Ce);
        try
            M{i}.iCe=M{i}.iCe;
        catch
            M{i}.iCe = inv(M{i}.Ce);
        end
    end
end


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
ipC = inv(M{2}.pC);
log_prior_t2 = spm_logdet(ipC)/2-0.5*Np*log(2*pi);

M{2}.ipC=ipC;
M{2}.log_prior_t2=log_prior_t2;
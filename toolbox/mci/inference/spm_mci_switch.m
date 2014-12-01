function [logp,logq1,logq2,logp1,logp2] = spm_mci_switch (P,M,U,Y,beta)
% Return log probability of tempered model switch
% FORMAT [logp,logq1,logq2,logp1,logp2] = spm_mci_switch (P,M,U,Y,beta)
%
% P,M,U,Y   as usual
% beta      inverse temperature (set to 1 to get usual posterior)
%
% logp      log tempered posterior probability
% logp1     log tempered posterior probability of model 1
% logp2     log tempered posterior probability of model 2
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_switch.m 6275 2014-12-01 08:41:18Z will $

loglike1= feval(M{1}.L,P,M{1},U{1},Y);
%p = M{1}.V'*(P-M{1}.pE);
p = P(:)-M{1}.V'*spm_vec(M{1}.pE);

logprior1=- p'*M{1}.ipC*p/2 + M{1}.log_prior_t2;

% Extract subset of parameters relevant to model 2
%P=P(M{2}.subset);
%p=M{1}.V'*(P-M{2}.pE);
p=P(:)- M{1}.V'*spm_vec(M{2}.pE);

% Subset of params for likelihood
P=P(M{2}.subset);
loglike2= feval(M{2}.L,P,M{2},U{2},Y);

% All params for prior
logprior2=- p'*M{2}.ipC*p/2 + M{2}.log_prior_t2;

% For sampling - must constrain extra params toward zero
logq1=loglike1+logprior1;
logq2=loglike2+logprior2;
logp=(1-beta)*logq1+beta*logq2;

% For scoring of models - ignore extra params
ind=M{2}.subset;
%p=(P-M{2}.pE(ind));
pE=spm_vec(M{2}.pE);
p=P-pE(ind);
logprior2=- p'*M{2}.ipC(ind,ind)*p/2 + M{2}.log_prior_t2;
logq2=loglike2+logprior2;

logp1=logprior1;
logp2=logprior2;

%logp1=loglike1;
%logp2=loglike2;
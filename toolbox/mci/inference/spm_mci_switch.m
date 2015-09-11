function [logp,logq1,logq2] = spm_mci_switch (Pr,M,U,Y,beta)
% Return log probability of tempered model switch
% FORMAT [logp,logq1,logq2] = spm_mci_switch (Pr,M,U,Y,beta)
%
% Pr        parameters (vectorised and in M.V subspace)
% M,U,Y     as usual
% beta      inverse temperature (set to 1 to get usual posterior)
%
% logp      log prob of model switch
% logq1     log joint of model 1
% logq2     log joint of model 2
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_switch.m 6548 2015-09-11 12:39:47Z will $

logq1 = spm_mci_joint(Pr,M{1},U{1},Y);

% Parameters in original space
P = M{1}.V*Pr+M{1}.vpE;
[L2,tmp,st] = feval(M{2}.L,P,M{2},U{2},Y);
e = P-M{2}.vpE;
L1 = - e'*M{2}.ipC*e/2 + M{2}.log_prior_t2;
logq2 = L1+L2;

logp=(1-beta)*logq1+beta*logq2;

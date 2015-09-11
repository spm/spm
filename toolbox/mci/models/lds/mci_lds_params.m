function [P] = mci_lds_params (M,U)
% LDS constrained: sample params from prior
% FORMAT [P] = mci_lds_params (M,U)
%
% M     Model structure
% U     Inputs
%
% P     Parameters
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_params.m 6548 2015-09-11 12:39:47Z will $

P=spm_normrnd(M.pE,M.pC,1);

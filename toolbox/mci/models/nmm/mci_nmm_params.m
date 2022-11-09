function [P] = mci_nmm_params (M,U)
% Generate parameters for two region NMM 
% FORMAT [P] = mci_nmm_params (M,U)
%
% M         Model structure
% U         Inputs
%
% P         Parameters
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

P=spm_normrnd(M.pE,M.pC,1);

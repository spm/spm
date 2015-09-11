function [P] = mci_nmm_params (M,U)
% Generate parameters for two region NMM 
% FORMAT [P] = mci_nmm_params (M,U)
%
% M         Model structure
% U         Inputs
%
% P         Parameters
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_nmm_params.m 6548 2015-09-11 12:39:47Z will $

P=spm_normrnd(M.pE,M.pC,1);

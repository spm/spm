function [Y] = mci_nmm_gen (M,U,P)
% Generate data from two region NMM 
% FORMAT [Y] = mci_nmm_gen (M,U,P)
%
% M         Model structure
% U         Inputs
% P         Parameters
%
% Y         Data
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

[G,x] = spm_mci_fwd (P,M,U);
e=randn(M.N,M.Nr)*sqrt(M.Ce);
Y=G+e;

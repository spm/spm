function [Y] = mci_ramsay_gen (P,M,U)
% Generate data from Ramsay model
% FORMAT [Y] = mci_ramsay_gen (P,M,U)
%
% P         Parameters
% M         Model structure
% U         Inputs
%
% Y         Data
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

[G,x] = spm_mci_fwd (P,M,U);
[N,l] = size(G);
e = randn(N,l)*sqrt(M.Ce);
Y = G + e;

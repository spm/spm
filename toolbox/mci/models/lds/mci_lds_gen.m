function [Y] = mci_lds_gen (M,U,P)
% LDS constrained: generate data
% FORMAT [Y] = mci_lds_gen (M,U,P)
%
% M     Model structure
% U     Inputs
% P     Parameters
%
% Y     Data
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

G = spm_mci_fwd (P,M,U);
e = M.sd*randn(M.N,M.d);
Y = G+e;


function [g,y] = mci_discount_gen (P,M,U)
% Output of discounting model
% FORMAT [g,y] = mci_discount_gen (P,M,U)
%
% P         parameters
% M         model structure
% U         U.X contains design matrix
%
% g         probability of taking option 1
% y         binary decisions based on g
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

a = mci_discount_act (P,M,U);
g = 1./(1+exp(-a));

y = g > rand(M.N,1);
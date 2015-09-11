function [g,y] = mci_logistic_gen (P,M,U)
% Output of logistic regression model
% FORMAT [g,y] = mci_logistic_gen (P,M,U)
%
% P         parameters
% M         model structure
% U         U.X contains design matrix
%
% g         probabilities of y=1
% y         binary decisions based on g
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_logistic_gen.m 6548 2015-09-11 12:39:47Z will $

a = mci_logistic_act (P,M,U);
g = 1./(1+exp(-a));

y = g > rand(M.N,1);
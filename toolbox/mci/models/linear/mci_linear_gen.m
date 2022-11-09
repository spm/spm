function [y] = mci_linear_gen (theta,M,U)
% Output of linear model
% FORMAT [y] = mci_linear_gen (theta,M,U)
%
% theta     regression coefficients
% M         model structure
% U         U.X contains design matrix
%
% y         outputs
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

y=U.X*theta(:);

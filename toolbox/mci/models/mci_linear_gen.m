function [y] = mci_linear_gen (theta,M,U)
% Output of linear model
% FORMAT [y] = mci_linear_gen (theta,M,U)
%
% theta     regression coefficients
% M         model structure
% U         U.X contains design matrix
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linear_gen.m 6275 2014-12-01 08:41:18Z will $

y=U.X*theta(:);

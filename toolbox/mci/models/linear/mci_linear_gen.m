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
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linear_gen.m 6548 2015-09-11 12:39:47Z will $

y=U.X*theta(:);

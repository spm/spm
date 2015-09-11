function [y] = mci_linsqr_gen (theta,M,U)
% Output of linear model with squared params
% FORMAT [y] = mci_linsqr_gen (theta,M,U)
%
% theta     regression coefficients
% M         model structure
% U         U.X contains design matrix
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

% Square parameters 
beta=theta(:).^2;

y=U.X*beta;


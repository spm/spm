function [L,E,st] = mci_linear_like (theta,M,U,Y)
% Compute log likelihood of linear model
% FORMAT [L,E,st] = mci_linear_like (theta,M,U,Y)
%
% theta     regression coefficients
% M         model
% U         inputs
% Y         data
% 
% L         Log likelihood
% E         Errors
% st        Status flag (0 for OK, -1 for problem)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linear_like.m 6275 2014-12-01 08:41:18Z will $

st=0;

yhat=U.X*theta(:);
T=length(yhat);
if isstruct(Y)
    E=sum(sum((Y.y-yhat).^2));
else
    E=sum(sum((Y-yhat).^2));
end

L = M.logdet_Ce - 0.5*T*log(2*pi);
L = L - 0.5*M.iCe*E;

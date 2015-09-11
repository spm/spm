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
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linear_like.m 6548 2015-09-11 12:39:47Z will $

st=0;

yhat=U.X*theta(:);
T=length(yhat);
if isstruct(Y)
    E=sum(sum((Y.y-yhat).^2));
else
    E=sum(sum((Y-yhat).^2));
end

L = -0.5*T*M.logdet_Ce - 0.5*T*log(2*pi);
L = L - 0.5*M.iCe*E;

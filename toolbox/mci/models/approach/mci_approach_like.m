function [L,yhat,st] = mci_approach_like (P,M,U,Y)
% Log-likelihood for approach model 
% FORMAT [L,yhat,st] = mci_approach_like (P,M,U,Y)
%
% P         parameters
% M,U,Y     as usual
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_approach_like.m 6548 2015-09-11 12:39:47Z will $

% Status flag (only used for dynamic systems)
st=[];

T=length(Y);
yhat = mci_approach_gen (P,M,U);
E=sum(sum((Y-yhat).^2));

L = -0.5*T*M.logdet_Ce - 0.5*T*log(2*pi);
L = L - 0.5*M.iCe*E;


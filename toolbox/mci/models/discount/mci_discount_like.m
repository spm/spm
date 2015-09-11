function [L,E,st] = mci_discount_like (P,M,U,Y)
% Compute log likelihood of discount model
% FORMAT [L,E,st] = mci_discount_like (P,M,U,Y)
%
% P         parameters
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
% $Id: mci_discount_like.m 6548 2015-09-11 12:39:47Z will $

st=0;

g = mci_discount_gen (P,M,U);

T=length(g);
if isstruct(Y)
    y=Y.y;
else
    y=Y;
end

L=0;
for n=1:T,
    L = L + y(n)*log(g(n)+eps)+(1-y(n))*log(1-g(n)+eps);
end
E=-L;




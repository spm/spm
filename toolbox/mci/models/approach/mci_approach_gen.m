function [y] = mci_approach_gen (P,M,U)
% Approach to limit model
% FORMAT [y] = mci_approach_gen (P,M,U)
%
% P         parameters
% M,U       as usual
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

V=exp(P(1));
tau=exp(P(2));
t=U.X;

y=-60+V*(1-exp(-t/tau));


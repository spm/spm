function [y] = mci_pb_gen (P,M,U)
% Preece-Baines growth model
% FORMAT [y] = mci_pb_gen (P,M,U)
%
% P         parameters
% M         model
% U         inputs
%
% y         time series
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_pb_gen.m 6548 2015-09-11 12:39:47Z will $

t=U.X;

s0=exp(P(1)); % Pre-pubertal velocity
s1=exp(P(2)); % Post-pubertal velocity
t0=P(3); % Age at growth spurt

h1=P(4); % asymptotic height
h0=P(5); % height at growth spurt

e1=exp(s0*(t-t0));
e2=exp(s1*(t-t0));

y=h1-2*(h1-h0)./(e1+e2);
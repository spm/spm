function [a] = mci_logistic_act (P,M,U)
% Activations of logistic model
% FORMAT [a] = mci_logistic_act (P,M,U)
%
% P         parameters
% M         model structure
% U         contains rewards and times
%
% a         activations of logistic model
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_logistic_act.m 6548 2015-09-11 12:39:47Z will $

a = U.X*P;
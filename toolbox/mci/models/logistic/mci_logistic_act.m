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

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

a = U.X*P;
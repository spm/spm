function [Y] = mci_rphase_gen (P,M,U)
% Generate data from reduced WCO model
% FORMAT [Y] = mci_rphase_gen (P,M,U)
%
% P     parameters
% M     model structure
% U     inputs
%
% Y     data
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

Y = spm_mci_fwd (P,M,U); 
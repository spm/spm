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
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

Y = spm_mci_fwd (P,M,U); 
function [dfdp] = mci_rphase_dfdp (x,u,P,M)
% Parameter sensitivity for phase model
% FORMAT [dfdp] = mci_rphase_dfdp (x,u,P,M)
%
% x      State vector
% u      inputs
% P      parameter vector
% M      model structure
%
% dfdp   Jacobian wrt. parameters, df/dp
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_rphase_dfdp.m 6548 2015-09-11 12:39:47Z will $

dfdp = spm_diff(M.f,x,u,P,M,3);



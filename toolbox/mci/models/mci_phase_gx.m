function [y,L] = mci_phase_gx (x,u,P,M)
% Observation function for phase model
% FORMAT [y,L] = mci_phase_gx (x,u,P,M)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_phase_gx.m 6275 2014-12-01 08:41:18Z will $

L=eye(M.n);
y=L*x;
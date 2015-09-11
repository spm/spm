function [y,L] = mci_phase_gx (x,u,P,M)
% Observation function for phase model
% FORMAT [y,L] = mci_phase_gx (x,u,P,M)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_phase_gx.m 6548 2015-09-11 12:39:47Z will $

L=eye(M.n);
y=L*x;
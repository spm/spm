function [y,L] = mci_phase_gx (x,u,P,M)
% Observation function for phase model
% FORMAT [y,L] = mci_phase_gx (x,u,P,M)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id$

L=eye(M.n);
y=L*x;
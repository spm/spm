function [y,L] = mci_phase_gx (x,u,P,M)
% Observation function for phase model
% FORMAT [y,L] = mci_phase_gx (x,u,P,M)
%__________________________________________________________________________

% Will Penny and Biswa Sengupta
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

L=eye(M.n);
y=L*x;
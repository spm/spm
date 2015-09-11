function [F] = mci_nmm_r2p2_dfdp (x,u,P,M)
% State Jacobian for two region, two parameter NMM
% FORMAT [F] = mci_nmm_r2p2_dfdp (x,u,P,M)
%
% x         State
% u         Inputs
% P         Parameters
% M         Model structure
%
% F         F(i,j) = df(x)_i/dtheta_j
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_nmm_r2p2_dfdx.m 6548 2015-09-11 12:39:47Z will $

% 18 state variables
F = zeros(18,18);

% f         Flow, dx/dt

curr_P=M.can_P; % Canonical parameter set

% 2 free parameters
curr_P.A{1}(2,1)=P(1); % Forward connection, w_21
curr_P.A{2}(1,2)=P(2); % Backward connection, w_12

f = mci_nmm_fx_delay(x,u,curr_P,M);


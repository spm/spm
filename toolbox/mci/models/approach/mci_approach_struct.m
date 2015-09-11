function [M,U] = mci_approach_struct (Nobs)
% Approach model structure
% FORMAT [M,U] = mci_approach_struct (Nobs)
%
% Nobs      Number of observations
% M         Model structure
% U         Input structure
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_approach_struct.m 6548 2015-09-11 12:39:47Z will $

M.l=1; % Single output variable

M.T=40;
dt=M.T/(Nobs-1);
U.X=[0:dt:M.T]';
M.N=length(U.X);
M.t=U.X;

M.L='mci_approach_like';
M.IS='mci_approach_gen';
M.dL='mci_approach_deriv';

M.pE=[log(20),log(5)]';
M.pC=[1/16 0; 0 1/16];

sigma_e=1;
M.Ce=sigma_e^2;
M.logdet_Ce=spm_logdet(M.Ce);
M.iCe=inv(M.Ce);




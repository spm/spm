function [M,U] = mci_pb_struct (Nobs)
% Preece-Baines model structure
% FORMAT [M,U] = mci_pb_struct (Nobs)
%
% Nobs      Number of observations
%
% M         Model structure
% U         Input structure
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_pb_struct.m 6548 2015-09-11 12:39:47Z will $

M.l=1; % Single output variable

M.T=30;
dt=M.T/(Nobs-1);
U.X=[0:dt:M.T]';
M.N=length(U.X);
M.t=U.X;

M.L='mci_pb_like';
M.IS='mci_pb_gen';
M.dL='mci_pb_deriv';

M.pE=[log(0.02),log(0.4),11,90,60]';
M.pC=diag(abs(M.pE/10));

sigma_e=1;
M.Ce=sigma_e^2;
M.logdet_Ce=spm_logdet(M.Ce);
M.iCe=inv(M.Ce);




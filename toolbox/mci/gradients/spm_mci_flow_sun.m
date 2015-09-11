function [f, flag, new_data] = spm_mci_flow_sun (t, x, data)
% Evaluate flow for Sundials routines
% FORMAT [f, flag, new_data] = spm_mci_flow_sun (t, x, data)
%
% t     time
% x     state
% data  .U inputs, .P parameters, .M model
%
% f     flow, dx/dt
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_flow_sun.m 6548 2015-09-11 12:39:47Z will $

U=data.U;
P=data.P;
M=data.M;

P=spm_unvec(P,M.pE);

f = spm_mci_flow_t(t,x,U,P,M);
flag = 0;
new_data = [];

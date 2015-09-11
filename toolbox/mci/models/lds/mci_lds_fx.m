function [f,A,Pt] = mci_lds_fx (x,u,P,M)
% Flow for linear system, dx/dt=Ax, with constrained connectivity
% FORMAT [f,A,Pt] = mci_lds_fx (x,u,P,M)
%
% x     State vector
% u     input
% P     parameters (vectorised)
% M     model structure
%
% f     Flow, dx/dt
% A     f=Ax
% Pt    Parameters (transformed from latent pars)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_fx.m 6548 2015-09-11 12:39:47Z will $

[A,Pt] = mci_lds_dfdx (x,u,P,M);
f=A*x;
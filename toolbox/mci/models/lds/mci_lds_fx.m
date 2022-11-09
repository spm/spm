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

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

[A,Pt] = mci_lds_dfdx (x,u,P,M);
f=A*x;
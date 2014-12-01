function [A,Pt] = mci_lds_dfdx (x,u,P,M)
% Jacobian for linear system, dx/dt=Ax, with constrained connectivity
% FORMAT [A,Pt] = mci_lds_dfdx (x,u,P,M)
%
% x     State vector
% u     input
% P     parameters (vectorised)
% M     model structure
%
% A     f=Ax
% Pt    Parameters (transformed from latent pars)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_dfdx.m 6275 2014-12-01 08:41:18Z will $

[Pt,a,b] = mci_lds_lat2par (P,M);
A=diag(a);

Nb=length(b);
for k=1:Nb,
    i=M.Aconn(k,1);
    j=M.Aconn(k,2);
    A(i,j)=b(k);
end


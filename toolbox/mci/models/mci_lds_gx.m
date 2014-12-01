function [y,L] = mci_lds_gx (x,u,P,M)
% Observation function for LDS
% FORMAT [y,L] = mci_lds_gx (x,u,P,M)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_gx.m 6275 2014-12-01 08:41:18Z will $

L=eye(M.n);
y=L*x;
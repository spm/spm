function [Pt,a,b] = mci_lds_lat2par (P,M)
% Convert latent params to params
% FORMAT [Pt,a,b] = mci_lds_lat2par (P,M)
%
% P         Parameters (latent)
% M         model structure
%
% Pt        Parameters (transformed)
% a         diagonal values
% b         off-diagonal values
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_lat2par.m 6548 2015-09-11 12:39:47Z will $

% Diagonal entries
s=exp(P(1:M.d));
a=s*M.a_typical;

% Off-diagonal entries
b=P(M.d+1:end);
    
Pt=[a(:);b(:)];
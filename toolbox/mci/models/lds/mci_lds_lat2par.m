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

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Diagonal entries
s=exp(P(1:M.d));
a=s*M.a_typical;

% Off-diagonal entries
b=P(M.d+1:end);
    
Pt=[a(:);b(:)];
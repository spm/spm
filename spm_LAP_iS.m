function [iS] = spm_LAP_iS(p,R)
% Default precision function for LAP models (hidden states)
% FORMAT [iS] = spm_LAP_iS(p,R)
%
% p{1} - vector of precisions for causal states (v)
% p{2} - vector of precisions for hidden states (v)
% R    - generalised precision matrix
%
% iS   - precision matrix for generalised states (causal and then hidden)
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging
 

% precisions for generalised errors on causal (V) and hidden states (W)
%==========================================================================
V   = kron(R,diag(exp(p.h)));
W   = kron(R,diag(exp(p.g)));
iS  = blkdiag(V,W);

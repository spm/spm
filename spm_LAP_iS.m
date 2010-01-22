function [iS] = spm_LAP_iS(p,R)
% default precision function for LAP models (hidden states)
% FORMAT [iS] = spm_LAP_iS(p,R)
%
% p{1} - vector of precisions for causal states (v)
% p{2} - vector of precisions for hidden states (v)
% R    - generalised precision matrix
%
% iS   - precision matrix for generalised states (causal and then hidden)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_LAP_iS.m 3694 2010-01-22 14:16:51Z karl $
 
% precisions for generalised errors on causal (V) and hidden states (W)
%==========================================================================
V   = kron(R,diag(exp(p.h)));
W   = kron(R,diag(exp(p.g)));
iS  = blkdiag(V,W);

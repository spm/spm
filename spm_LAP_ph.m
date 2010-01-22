function p = spm_LAP_ph(x,v,h,M);
% default precision function for LAP models (causal states)
% FORMAT p = spm_LAP_ph(x,v,h,M);
%
% x  - hidden states
% v  - causal states
% h  - precision parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_LAP_ph.m 3694 2010-01-22 14:16:51Z karl $

% fixed components
%--------------------------------------------------------------------------
p = sparse(M.l,1);
try
    V = diag(M.V);
    if all(V)
        p = log(V);
    end
end

% free components
%--------------------------------------------------------------------------
for i = 1:length(M.Q)
    p = p + h(i)*diag(M.Q{i});
end

function p = spm_LAP_ph(x,v,h,M)
% Default precision function for LAP models (causal states)
% FORMAT p = spm_LAP_ph(x,v,h,M)
%
% x  - hidden states
% v  - causal states
% h  - precision parameters
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


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

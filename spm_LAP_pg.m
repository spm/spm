function p = spm_LAP_pg(x,v,h,M)
% Default precision function for LAP models (hidden states)
% FORMAT p = spm_LAP_pg(x,v,h,M)
%
% x  - hidden states
% v  - causal states
% h  - precision parameters
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


% fixed components
%--------------------------------------------------------------------------
p = sparse(M.n,1);
try
    W = diag(M.W);
    if all(W)
        p = log(W);
    end
end

% free components
%--------------------------------------------------------------------------
for i = 1:length(M.R)
    p = p + h(i)*diag(M.R{i});
end

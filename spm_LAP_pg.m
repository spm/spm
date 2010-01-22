function p = spm_LAP_pg(x,v,h,M);
% default precision function for LAP models (hidden states)
% FORMAT p = spm_LAP_pg(x,v,h,M);
%
% x  - hidden states
% v  - causal states
% h  - precision parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_LAP_pg.m 3694 2010-01-22 14:16:51Z karl $

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

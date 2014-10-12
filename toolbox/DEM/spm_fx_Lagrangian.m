function [f] = spm_fx_Lagrangian(P,M,U)
% FORMAT [f] = spm_fx_Lagrangian(P,M,U)
%
% flow subfunction for Langrangian demo



% Lagrangian and Q
%==========================================================================
[n m] = size(U);

for i = 1:m
    x = U(:,i);
    L = P.P{1};
    for j = 2:length(P.P)
        L = L + P.P{j}(:)'*x;
        x = kron(x,U(:,i));
    end
    
end

% estimate parameters of Lagrangian and Q
%==========================================================================
f    = P.Q*dLdx;




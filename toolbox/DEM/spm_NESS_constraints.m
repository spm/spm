function [ks,kq,kg,kh] = spm_NESS_constraints(o,A,K,L)
% constraints on polynomial coefficients or dynamical systems
% FORMAT [ks,kq,kg,kh] = spm_NESS_constraints(o,A,K,L);
% o - matrix of orders for polynomial expansion
% A - adjacency matrix (dynamical coupling)
% K - upper bound on order for surprisal parameters
% J - upper bound on order for flow operator parameters
%
% ks  - indices for surprisal   parameters
% kq  - indices for solenoidal  parameters
% kg  - indices for dissipative parameters
% kh  - indices for curvature   parameters

%
%--------------------------------------------------------------------------
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% constraints on potential parameters due to dissipative flow
%--------------------------------------------------------------------------
[n,nb] = size(o);                            % number of basis functions

ks    = sum(o) > K;                          % polynomial order constraints
ks    = ks | ~sum(o);                        % suppress constant
for i = 1:n
    for j = 1:n
        if ~A(i,j)
            ks = ks | (o(i,:) & o(j,:));
        end
    end
end

% diagonal terms of Hessian S
%--------------------------------------------------------------------------
kh    = any(o == 2);

% constraints on the order of the polynomial expansion
%--------------------------------------------------------------------------
k     = cell(n,n);
for i = 1:n
    for j = 1:n
        k{i,j} = sum(o) > L;
    end
end

% constraints due to diagonal elements of Hessian
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        if ~A(i,j)
            k{i,j} = ones(1,nb);
            k{j,i} = ones(1,nb);
        end
    end
end

% constraints due to non-negative gradients
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        if ~A(i,j)
            for q = 1:n
                k{i,q} = k{i,q} | o(j,:);
            end
        end
    end
end

% assemble and combine constraints for Q
%--------------------------------------------------------------------------
kq    = [];
for i = 1:n
    for j = i:n
        kq = [kq (k{i,j} | k{j,i})];
    end
end

% constraints due to diagonal elements of Hessian S
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        if i == j
            k{i,j} = ones(1,nb);
        else
            k{j,i} = zeros(1,nb);
        end
    end
end

kg    = [];
for i = 1:n
    for j = i:n
        kg = [kg (k{i,j} | k{j,i})];
    end
end


return






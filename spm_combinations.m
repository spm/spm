function U = spm_combinations(Nu)
% FORMAT U = spm_combinations(Nu)
% Nu  - vector of dimensions
% U   - combinations of indices
%
% returns a matrix of all combinations of Nu
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


Nf    = numel(Nu);
U     = zeros(prod(Nu),Nf);
for f = 1:Nf
    for j = 1:Nf
        if j == f
            k{j} = 1:Nu(j);
        else
            k{j} = ones(1,Nu(j));
        end
    end
    u     = 1;
    for i = 1:Nf
        u = kron(k{i},u);
    end

    % accumulate
    %----------------------------------------------------------------------
    U(:,f) = u(:);
    
end

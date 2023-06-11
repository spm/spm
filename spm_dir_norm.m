function A = spm_dir_norm(A)
% Normalisation of a (Dirichlet) conditional probability matrix
% FORMAT A = spm_dir_norm(a)
%
% a    - (Dirichlet) parameters of a conditional probability matrix
%
% A    - conditional probability matrix
%__________________________________________________________________________

% Karl Friston 
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% deal with cells
%--------------------------------------------------------------------------
if iscell(A)
    for g = 1:numel(A)
        A{g} = spm_dir_norm(A{g});
    end
    return
end

% deal with Dirichlet tensors
%--------------------------------------------------------------------------
A0      = sum(A,1);
i       = logical(A0);
A       = rdivide(A,A0);
A(:,~i) = 1/size(A,1);
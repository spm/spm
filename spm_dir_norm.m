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

A           = rdivide(A,sum(A,1));
A(isnan(A)) = 1/size(A,1);

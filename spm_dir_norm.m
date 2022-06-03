function A  = spm_dir_norm(A)
% normalisation of a (Dirichlet) conditional probability matrix
% FORMAT A  = spm_dir_norm(a)
%
% a    - (Dirichlet) parameters of a conditional probability matrix
%
% A    - conditional probability matrix
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston 
% $Id: spm_KL_dir.m 8183 2021-11-04 15:25:19Z guillaume $

%--------------------------------------------------------------------------
A           = bsxfun(@rdivide,A,sum(A,1));
A(isnan(A)) = 1/size(A,1);
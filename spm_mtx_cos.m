function c = spm_mtx_cos(A,B)
% returns the cosine of the angle between A and B
% FORMAT c = spm_mtx_cos(A,B)
%
% a    - (Dirichlet) parameters of a conditional probability matrix
%
% c = arccos( <A|B> /(<A|A><B|B>))
%__________________________________________________________________________

% Karl Friston 
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

A = A(:);
B = B(:);
c = acos( A'*B / (sqrt(A'*A)*sqrt(B'*B)));
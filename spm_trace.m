function [C] = spm_trace(A,B)
% Fast trace for large matrices: C = spm_trace(A,B) = trace(A*B)
% FORMAT [C] = spm_trace(A,B)
%
% C = spm_trace(A,B) = trace(A*B) = sum(sum(A'.*B));
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging

% fast trace for large matrices: C = spm_trace(A,B) = trace(A*B)
%--------------------------------------------------------------------------
C = sum(sum(A'.*B));

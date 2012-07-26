function [C] = spm_trace(A,B)
% fast trace for large matrices: C = spm_trace(A,B) = trace(A*B)
% FORMAT [C] = spm_trace(A,B)
%
% C = spm_trace(A,B) = trace(A*B) = sum(sum(A'.*B));
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_trace.m 4805 2012-07-26 13:16:18Z karl $

% fast trace for large matrices: C = spm_trace(A,B) = trace(A*B)
%--------------------------------------------------------------------------
C = sum(sum(A'.*B));

function ndx = spm_sub2ind(siz,v1,v2,varargin)
%SUB2IND Linear index from multiple subscripts (without error checking).
%   SUB2IND is used to determine the equivalent single index
%   corresponding to a given set of subscript values.
%
%   IND = SUB2IND(SIZ,I,J) returns the linear index equivalent to the
%   row and column subscripts in the arrays I and J for a matrix of
%   size SIZ. 
%
%   IND = SUB2IND(SIZ,I1,I2,...,IN) returns the linear index
%   equivalent to the N subscripts in the arrays I1,I2,...,IN for an
%   array of size SIZ.
%
%   I1,I2,...,IN must have the same size, and IND will have the same size
%   as I1,I2,...,IN. For an array A, if IND = SUB2IND(SIZE(A),I1,...,IN)),
%   then A(IND(k))=A(I1(k),...,IN(k)) for all k.
%
%   Class support for inputs I,J: 
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also IND2SUB.

%   Copyright 1984-2015 The MathWorks, Inc.

siz    = double(siz);
lensiz = length(siz);

numOfIndInput = nargin - 1;
if lensiz < numOfIndInput
    % Adjust for trailing singleton dimensions
    siz = [siz, ones(1,numOfIndInput - lensiz)];
elseif lensiz > numOfIndInput
    % Adjust for linear indexing on last element
    siz = [siz(1:numOfIndInput-1), prod(siz(numOfIndInput:end))];
end

% Compute linear indices
%--------------------------------------------------------------------------
ndx   = double(v1);
if numOfIndInput >= 2
    ndx = ndx + (double(v2) - 1).*siz(1);
end 
if numOfIndInput > 2
    k = cumprod(siz);
    for i = 3:numOfIndInput
        v   = varargin{i-2};
        ndx = ndx + (double(v) - 1)*k(i - 1);
    end
end
function t = spm_convmtx(v,n)
% as for convmtx but for sparse matrices
% FORMAT t = spm_convmtx(v,n)
%
%--------------------------------------------------------------------------
%   CONVMTX(C,N) returns the convolution matrix for vector C.
%   If C is a column vector and X is a column vector of length N,
%   then CONVMTX(C,N)*X is the same as CONV(C,X).
%   If R is a row vector and X is a row vector of length N,
%   then X*CONVMTX(R,N) is the same as CONV(R,X).
%   See also CONV.
%_______________________________________________________________________
% Copyright (C) 1988-2004 The MathWorks, Inc.

% L. Shure and T. Krauss
% $Id: spm_convmtx.m 3589 2009-11-20 17:17:41Z guillaume $

[mv,nv] = size(v);
v = v(:);                                     % make v a column vector
c = [v; sparse(n-1,1)];
r = sparse(n,1);
m = length(c);
x = [r(n:-1:2) ; c(:)];                       % build vector of user data
cidx = (0:m-1)';
ridx = n:-1:1;
t = cidx(:,ones(n,1)) + ridx(ones(m,1),:);    % Toeplitz subscripts
t(:) = x(t);                                  % actual data
% end of toeplitz code

if mv < nv
    t = t.';
end


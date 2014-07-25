function t = spm_convmtx(v,n,OPT)
% as for convmtx but with boundary conditions
% FORMAT t = spm_convmtx(C,N,OPT)
%
% OPT  - 'circular' boundary conditions
%      - 'square'   top and tail convolution matrix
%
%--------------------------------------------------------------------------
%   CONVMTX(C,N) returns the convolution matrix for vector C.
%   If C is a column vector and X is a column vector of length N,
%   then CONVMTX(C,N)*X is the same as CONV(C,X).
%   If R is a row vector and X is a row vector of length N,
%   then X*CONVMTX(R,N) is the same as CONV(R,X).
%   See also CONV.%
%   With the circular option the convolution matrix is reduced to N X N
%__________________________________________________________________________
% Copyright (C) 1988-2004 The MathWorks, Inc.
 
% L. Shure and T. Krauss
% $Id: spm_convmtx.m 6122 2014-07-25 13:48:47Z karl $
 
if nargin < 3;
    OPT = 'none';
end
 
% create Toeplitz matrix
%--------------------------------------------------------------------------
[mv,nv] = size(v);
v    = full(v(:));                              % make v a column vector
c    = [v; zeros(n-1,1)];
r    = zeros(n,1);
m    = length(c);
x    = [r(n:-1:2) ; c(:)];                      % build vector of user data
cidx = (0:m-1)';
ridx = n:-1:1;
t    = cidx(:,ones(n,1)) + ridx(ones(m,1),:);   % Toeplitz subscripts
t(:) = x(t);                                    % actual data
 
% transpose if necessary
%--------------------------------------------------------------------------
if mv < nv
    t = t.';
end
 
% apply optional boundary conditions
%--------------------------------------------------------------------------
switch OPT
    
    case('circular')
        m      = fix((size(t,1) - n)/2);
        j      = (1:m) + m;
        t(j,:) = t(j,:) + t(j + n,:);
        j      = (1:m) + n;
        t(j,:) = t(j,:) + t(j - n,:);
        j      = 1:n;
        t      = t(j + m,:);
        
    case('square')
        m      = fix((size(t,1) - n)/2);
        j      = 1:n;
        t      = t(j + m,:);
        
    otherwise
end

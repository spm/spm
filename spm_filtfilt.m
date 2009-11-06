function y = spm_filtfilt(b,a,x)
% Zero-phase forward and reverse digital filtering
% FORMAT y = spm_filtfilt(b,a,x)
%
% b        - filter parameters (numerator)
% a        - filter parameters (denominator)
% x        - input data vector (if matrix, filter over columns)
%
% y        - filtered data
%__________________________________________________________________________
%
% The filter is described by the difference equation:
%
%     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
% After filtering in the forward direction, the filtered sequence is then
% reversed and run back through the filter; Y is the time reverse of the
% output of the second filtering operation. The result has precisely zero
% phase distortion and magnitude modified by the square of the filter's
% magnitude response.  Care is taken to minimize startup and ending
% transients by matching initial conditions.
%
% The length of the input x must be more than three times
% the filter order, defined as max(length(b)-1,length(a)-1).
%
% References:
% [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed, McGraw-Hill, 2001
% [2] Fredrik Gustafsson, Determining the initial states in forward-backward
%     filtering, IEEE Transactions on Signal Processing, pp. 988--992,
%     April 1996, Volume 44, Issue 4
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
 
% L. Shure, T. Krauss, F. Gustafsson
% Copyright 1988-2004 The MathWorks, Inc.
% $Id: spm_filtfilt.m 3541 2009-11-06 17:34:40Z guillaume $

% Check input data
%--------------------------------------------------------------------------
[m,n] = size(x);
if n>1 && m>1
    y = zeros(size(x));
    for i=1:n
        y(:,i) = spm_filtfilt(b,a,x(:,i));
    end
    return
end
if m==1, x = x(:); end
len = size(x,1);

% Check filter parameters
%--------------------------------------------------------------------------
b  = b(:).';    a  = a(:).';
nb = length(b); na = length(a);
nfilt = max(nb,na);
if nb < nfilt, b(nfilt)=0; end
if na < nfilt, a(nfilt)=0; end

nfact = 3*(nfilt-1);
if len <= nfact
    error('Data must have length more than 3 times filter order.');
end

% Use sparse matrix to solve system of linear equations for initial 
% conditions zi are the steady-state states of the filter b(z)/a(z) in the
% state-space implementation of the 'filter' command
%--------------------------------------------------------------------------
rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
data = [1+a(2) a(3:nfilt) ones(1,nfilt-2) -ones(1,nfilt-2)];
sp   = sparse(rows,cols,data);
zi   = sp \ (b(2:nfilt).' - a(2:nfilt).'*b(1));

% Extrapolate beginning and end of data sequence using a "reflection
% method".  Slopes of original and extrapolated sequences match at the end
% points. This reduces end effects
%--------------------------------------------------------------------------
y    = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];

% Filter, reverse data, filter again, and reverse data again
%--------------------------------------------------------------------------
y    = filter(b,a,y,zi*y(1));
y    = y(end:-1:1);
y    = filter(b,a,y,zi*y(1));
y    = y(end:-1:1);

% Reformat y
%--------------------------------------------------------------------------
y([1:nfact len+nfact+(1:nfact)]) = [];
if m == 1, y = y.'; end

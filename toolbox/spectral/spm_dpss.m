function [E] = spm_dpss (N,NW)
% Compute discrete prolate spheroidal sequences
% FORMAT [E] = spm_dpss (N,NW)
%
% N         Length of taper
% NW        Product of N and W
% 
% E         [N x 2NW] matrix containing dpss sequences
%           The kth column contains the sequence which
%           comprises the length N signal that is kth most 
%           concentrated in the frequency band |w|<=2*pi*W radians 
%
% See Section 8.3 in
%     Percival, D.B. and Walden, A.T., "Spectral Analysis For Physical
%     Applications", Cambridge University Press, 1993. 
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Eric Breitenberger, 10/3/95
% Copyright 1988-2004 The MathWorks, Inc.
% $Id: spm_dpss.m 3821 2010-04-15 14:26:34Z will $

W=NW/N;

k = min(round(2*N*W),N);
k = max(k,1);
k = [1 k];

% Generate first and second diagonals
d=((N-1-2*(0:N-1)').^2)*.25*cos(2*pi*W); 
ee=(1:N-1)'.*(N-1:-1:1)'/2;              

% Get the eigenvalues of B - see page 382 Percival and Walden
v = tridieig(d,[0; ee],N-k(2)+1,N-k(1)+1);
v = v(end:-1:1);
Lv = length(v);

% Compute the eigenvectors by inverse iteration with
% starting vectors of roughly the right shape.
E = zeros(N,k(2)-k(1)+1);
t = (0:N-1)'/(N-1)*pi;
for j = 1:Lv
   e = sin((j+k(1)-1)*t);
   
   e = tridisolve(ee,d-v(j),ee,e);
   e = tridisolve(ee,d-v(j),ee,e/norm(e));
   e = tridisolve(ee,d-v(j),ee,e/norm(e));

   E(:,j) = e/norm(e);
end

d=mean(E);
for i=k(1):k(2)
   if rem(i,2)  % i is odd
     % Polarize symmetric dpss
       if d(i-k(1)+1)<0, E(:,i-k(1)+1)=-E(:,i-k(1)+1); end
   else         % i is even
     % Polarize anti-symmetric dpss
       if E(2,i-k(1)+1)<0, E(:,i-k(1)+1)=-E(:,i-k(1)+1); end
   end
end


function x = tridieig(c,b,m1,m2,eps1);
%TRIDIEIG  Find a few eigenvalues of a tridiagonal matrix.
%   LAMBDA = TRIDIEIG(D,E,M1,M2).  D and E, two vectors of length N,
%   define a symmetric, tridiagonal matrix:
%      A = diag(E(2:N),-1) + diag(D,0) + diag(E(2:N),1)
%   E(1) is ignored.
%   TRIDIEIG(D,E,M1,M2) computes the eigenvalues of A with indices
%      M1 <= K <= M2.
%   TRIDIEIG(D,E,M1,M2,TOL) uses TOL as a tolerance.
%
%   C. Moler
%   Copyright 1988-2002 The MathWorks, Inc.

if nargin < 5, eps1 = 0; end
n = length(c);
b(1) = 0;
beta = b.*b;
xmin = min(c(n) - abs(b(n)),min(c(1:n-1) - abs(b(1:n-1)) - abs(b(2:n))));
xmax = max(c(n) + abs(b(n)),max(c(1:n-1) + abs(b(1:n-1)) + abs(b(2:n))));
eps2 = eps*max(xmax,-xmin);
if eps1 <= 0, eps1 = eps2; end
eps2 = 0.5*eps1 + 7*eps2;

x0 = xmax;
x = zeros(n,1);
wu = zeros(n,1);
x(m1:m2) = xmax(ones(m2-m1+1,1));
wu(m1:m2) = xmin(ones(m2-m1+1,1));
z = 0;
for k = m2:-1:m1
   xu = xmin;
   for i = k:-1:m1
      if xu < wu(i)
         xu = wu(i);
         break
      end
   end
   if x0 > x(k), x0 = x(k); end
   while 1
      x1 = (xu + x0)/2;
      if x0 - xu <= 2*eps*(abs(xu)+abs(x0)) + eps1
         break
      end
      z = z + 1;
      a = 0;
      q = 1;
      for i = 1:n
         if q ~= 0
            s = beta(i)/q;
         else
            s = abs(b(i))/eps;
         end
         q = c(i) - x1 - s;
         a = a + (q < 0);
      end
      if a < k
         if a < m1
            xu = x1;
            wu(m1) = x1;
         else
            xu = x1;
            wu(a+1) = x1;
            if x(a) > x1, x(a) = x1; end
         end
      else
         x0 = x1;
      end
   end
   x(k) = (x0 + xu)/2;
end
x = x(m1:m2)';

function x = tridisolve(a,b,c,d)
%   TRIDISOLVE  Solve tridiagonal system of equations.
%     x = TRIDISOLVE(a,b,c,d) solves the system of linear equations
%     b(1)*x(1) + c(1)*x(2) = d(1),
%     a(j-1)*x(j-1) + b(j)*x(j) + c(j)*x(j+1) = d(j), j = 2:n-1,
%     a(n-1)*x(n-1) + b(n)*x(n) = d(n).
%
%   The algorithm does not use pivoting, so the results might
%   be inaccurate if abs(b) is much smaller than abs(a)+abs(c).
%   More robust, but slower, alternatives with pivoting are:
%     x = T\d where T = diag(a,-1) + diag(b,0) + diag(c,1)
%     x = S\d where S = spdiags([[a; 0] b [0; c]],[-1 0 1],n,n)

x = d;
n = length(x);

for j = 1:n-1
   mu = a(j)/b(j);
   b(j+1) = b(j+1) - mu*c(j);
   x(j+1) = x(j+1) - mu*x(j);
end

x(n) = x(n)/b(n);
for j = n-1:-1:1
   x(j) = (x(j)-c(j)*x(j+1))/b(j);
end



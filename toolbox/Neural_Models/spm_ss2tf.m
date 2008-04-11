function [num, den] = ss2tf(a,b,c,d,iu)
%SS2TF  State-space to transfer function conversion.
%   [NUM,DEN] = SS2TF(A,B,C,D,iu)  calculates the transfer function:
%
%               NUM(s)          -1
%       H(s) = -------- = C(sI-A) B + D
%               DEN(s)
%   of the system:
%       .
%       x = Ax + Bu
%       y = Cx + Du
%
%   from the iu'th input.  Vector DEN contains the coefficients of the
%   denominator in descending powers of s.  The numerator coefficients
%   are returned in matrix NUM with as many rows as there are 
%   outputs y.
%
%   See also TF2SS, ZP2TF, ZP2SS.

%   J.N. Little 4-21-85
%   Revised 7-25-90 Clay M. Thompson, 10-11-90 A.Grace
%   Copyright 1984-2003 The MathWorks, Inc. 
%   $Revision: 1.24.4.1 $  $Date: 2003/05/01 20:42:46 $

error(nargchk(4,5,nargin));
[msg,a,b,c,d]=abcdchk(a,b,c,d); error(msg);

[mc,nu] = size(d);
if nargin==4,
  if (nu<=1)
    iu = 1;
  else
    error('MATLAB:ss2tf:NeedIU',...
          'IU must be specified for systems with more than one input.');
  end
end

den = poly(a);
if ~isempty(b), b = b(:,iu); end
if ~isempty(d), d = d(:,iu); end

% System is just a gain or it has only a denominator:
if isempty(b) && isempty(c)
    num = d;
    if isempty(d) && isempty(a)
        den = [];
    end
    return;
end

nc = length(a);
num = ones(mc, nc+1);
for i=1:mc
    num(i,:) = poly(a-b*c(i,:)) + (d(i) - 1) * den;
end

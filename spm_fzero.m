function b = spm_fzero(FunFcn,x,tol,trace,p1,p2,p3,p4,p5,p6,p7,p8,p9)
% Find a zero of a function of one variable
% FORMAT b = spm_fzero(FunFcn,x,tol,trace,p1,p2,p3,p4,p5,p6,p7,p8,p9)
% FunFcn - String containing name of function
% x      - initial guess
% tol    - (optional) relative tolerance for the convergence test
%          (leave empty for default value)
% trace  - (optional) non-zero value triggers a printing trace of the steps
%          (leave empty for default value)
% p1-p9  - additional constant parameters, passed to function FunFcn.
%          spm_fzero finds zero of eval([FunFcn '(x,p1,p2,..)']) over
%          values of x.
%__________________________________________________________________________
%
% spm_fzero finds a zero of a function of one variable. The value
% returned is near a point where F changes sign. For example,
% spm_fzero('sin',3) is pi.  Note the quotes around sin. Ordinarily,
% functions are defined in M-files.
%
% Adapted from MathWorks function FZERO, which does not allow additional
% parameters to be passed to the function.
%
% This algorithm was originated by T. Dekker.  An Algol 60 version,
% with some improvements, is given by Richard Brent in "Algorithms for
% Minimization Without Derivatives", Prentice-Hall, 1973.  A Fortran
% version is in Forsythe, Malcolm and Moler, "Computer Methods
% for Mathematical Computations", Prentice-Hall, 1976.
%
%__________________________________________________________________________
% %W% Andrew Holmes, Loren Shure, Clive Moler %E%

% Version history
% - Clive Moler - 1/86
% - Clive Moler - 3/87
% - Loren Shure - 12/88
% - Andrew Holmes - 12/93

% Initialization

if nargin < 4, trace = []; end
if nargin < 3, tol = []; end
if nargin<2 error('insufficient arguments'), end
if isempty(trace) trace=0; end
if isempty(tol) tol=eps; end

% set up function call
args = [];
for n = 1:nargin-4, args = [args,',p',int2str(n)]; end
argsa = ['(a',args,')'];
argsb = ['(b',args,')'];

if trace, clc, end
if (length(x) > 1) | (~finite(x))
   error('Second argument must be a finite scalar.')
end
if x ~= 0, dx = x/20;
else, dx = 1/20;
end
a = x - dx;  fa = eval([FunFcn,argsa]);
if trace, home, init = [a fa], end
b = x + dx;  fb = eval([FunFcn,argsb]);
if trace, home, init = [b fb], end

% Find change of sign.

while (fa > 0) == (fb > 0)
   dx = 2*dx;
   a = x - dx;  fa = eval([FunFcn,argsa]);
   if trace, home, sign = [a fa], end
   if (fa > 0) ~= (fb > 0), break, end
   b = x + dx;  fb = eval([FunFcn,argsb]);
   if trace, home, sign = [b fb], end
end

fc = fb;
% Main loop, exit from middle of the loop
while fb ~= 0
   % Insure that b is the best result so far, a is the previous
   % value of b, and c is on the opposite of the zero from b.
   if (fb > 0) == (fc > 0)
      c = a;  fc = fa;
      d = b - a;  e = d;
   end
   if abs(fc) < abs(fb)
      a = b;    b = c;    c = a;
      fa = fb;  fb = fc;  fc = fa;
   end

   % Convergence test and possible exit
   m = 0.5*(c - b);
   toler = 2.0*tol*max(abs(b),1.0);
   if (abs(m) <= toler) + (fb == 0.0), break, end

   % Choose bisection or interpolation
   if (abs(e) < toler) + (abs(fa) <= abs(fb))
   % Bisection
      d = m;  e = m;
   else
   % Interpolation
      s = fb/fa;
      if (a == c)
      % Linear interpolation
         p = 2.0*m*s;
         q = 1.0 - s;
      else
      % Inverse quadratic interpolation
         q = fa/fc;
         r = fb/fc;
         p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
         q = (q - 1.0)*(r - 1.0)*(s - 1.0);
      end;
      if p > 0, q = -q; else p = -p; end;
      % Is interpolated point acceptable
      if (2.0*p < 3.0*m*q - abs(toler*q)) * (p < abs(0.5*e*q))
         e = d;  d = p/q;
      else
         d = m;  e = m;
      end;
   end % Interpolation

   % Next point
   a = b;
   fa = fb;
   if abs(d) > toler, b = b + d;
   else if b > c, b = b - toler;
        else b = b + toler;
        end
   end
   fb = eval([FunFcn,argsb]);
   if trace, home, step = [b fb], end
end % Main loop

function [x,a,b,i] = spm_fzeroIB(f,a,b,z,tol,varargin)
% Find a zero of a function of one variable using interval bisection
% FORMAT [x,a,b,i] = spm_fzeroIB(f,a,b,z,tol,P1,P2,...)
%
% f      - string naming real-valued function of single real variable
%          f must be vectorised, and may be an inline function
% a      - lower bounds for x s.t. f(x)=z (scalar or same size as z)
% b      - upper bounds for x s.t. f(x)=z (scalar or same size as z)
% z      - zero offset(s) - seek x('s) at which f(x)=z rather than f(x)=0
%          vector of offsets useful for finding multiple crossing points
% tol    - Tolerance [Default eps*10]
% P1,... - Optional dditional parameters passed to function, f(x,P1,P2,...)
%        - Pass empty matrix for tol to use default value.
% i      - Number of iterations
%_______________________________________________________________________
%
% spm_fzeroIB tries to find a zero of the real valued function f of a
% single variable, using interval bisection.
% 
% f can be a mFunction or an inline function, which must be
% vectorised.
% 
% Zero offset(s) z may be supplied, whence the value(s) sought is(are)
% x such that f(x) = z, i.e. f(x)-z=0. The algorithm is vectorised, so
% matrix of zero offsets may be given.
% 
% Starting estimates must be provided as an interval [a,b] within which
% f changes sign (or crosses z in the case of zero offsets). Thus
% f(a)-z & f(b)-z must differ in sign. For vector offsets z, a & b can
% be scalars, or they can be matrices of bounds corresponding to the
% offsets in z.
%
% a and b returned specify the interval(s), smaller than tol, which
% contains the point where f(x)=z. x is the midpoint of this interval.
%
% Additional parameters to be passed to f can be specified as
% P1,P2,..., such that the values of x where f(x,P1,P2,...)=z are
% sought.
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%-Parameters
%-----------------------------------------------------------------------
Dtol  = 10^-10;
maxIt = 10000;

%-Check arguments
%-----------------------------------------------------------------------
if (nargin<5 | isempty(tol)), tol = Dtol; end
if (nargin<4 | isempty(z)  ),   z = 0;    end
if nargin<3, error('Insufficient arguments'); end
if ~isfinite([a(:);b(:);z(:)]), error('ranges (& offsets) must be finite'), end

%-Computation - Interval Bisection
%=======================================================================
a = ones(size(z)) .* a;
b = ones(size(z)) .* b;
fa = feval(f,a,varargin{:})-z;
fb = feval(f,b,varargin{:})-z;

%-Check limits bound a zero
if any(fa(:).*fb(:)>0), error('Starting interval(s) don''t contain fzero'), end

%-Initialise
i = 0;
x = a+(b-a)/2;
Q = 1:prod(size(z));

%-Interval bisection
%-----------------------------------------------------------------------
while length(Q) &  i<maxIt
	fxQ       = feval(f,x(Q),varargin{:})-z(Q);
	mQ        = fa(Q).*fxQ > 0;
	a(Q(mQ))  = x(Q(mQ));   fa(Q(mQ))  = fxQ(mQ);
	b(Q(~mQ)) = x(Q(~mQ));	fb(Q(~mQ)) = fxQ(~mQ);
	x(Q)      = a(Q) + (b(Q)-a(Q))/2;
	Q         = Q( (b(Q)-a(Q))>tol );
	i         = i+1;
end

if i==maxIt, warning('convergence criteria not reached - maxIt reached'), end

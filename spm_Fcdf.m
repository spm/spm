function p = spm_Fcdf(f,df,OffSet)
% Cumulative Distribution Function (CDF) of F distribution
% FORMAT p = spm_Fcdf(f,df,[OffSet])
% f      - F-variate
% df     - (vector of) degrees of freedom
% OffSet - (optional) offset, subtracted from resulting probability.
%          This enables use with fzero to compute inverse CDF.
%
%        - Multiple f may be supplied, with either constant df, or matrix of
%          pairs of df corresponding to the elements of f(:)'
%        - df must be a 2 x n (or n x 2) matrix of pairs of degrees of
%          freedom. (2x2 matrices are read as having df pairs in columns)
%__________________________________________________________________________
%
% spm_Fcdf implements the Cumulative Distribution of the F-distributions.
%
% Returns the probability p, that an F-variate on df degrees of freedom
% is less than f. p = Pr(F <= f) for F, an F random variable on df degrees of 
% freedom.
%
% Uses the incomplete beta function, as described in "Numerical Recipies in C"
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version control-%
% V1b	- 24/08/94 - streamlined code

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<3 OffSet=0; end
if nargin<2 error('insufficient arguments'), end

if any(abs(OffSet(:)-0.5)>0.5) error('OffSet must be in [0,1]'), end
if any(df(:)<=0) error('df must be strictly positive'), end

%-re-orient df to 2 x n size, if n x 2 with n ~= 2.
%---------------------------------------------------------------------------
if (size(df,1) ~= 2) & (size(df,2) == 2) df = df'; end
if  size(df,1) ~= 2 error('df must have 2 rows or 2 columns'), end

%-check sizes of arguments
%---------------------------------------------------------------------------
if prod(size(f)) == size(df,2)
elseif length(f) == 1
	f  = f*ones(1,size(df,2));
elseif size(df,2) == 1
	df = meshgrid(df,1:prod(size(f)))';
else
	error('f and df not of compatible size'), end
end % if (size)

%-Store sizes-reshape arguments to rows
%---------------------------------------------------------------------------
[p_rdim,p_cdim] = size(f);
f   = f(:)';
n   = length(f);

%-Computation, avoiding cases f<=0 where p=0
%---------------------------------------------------------------------------
p=zeros(1,n);
K  = find(f>0);
if length(K)>0;
	df1K = df(1,K); df2K = df(2,K);
	p(K) = 1 - betainc(df2K./(df2K + df1K.*f(K)),df2K/2,df1K/2);
end % f>0 cases

p  = reshape(p,p_rdim,p_cdim); end

%-Subtract (which defaults to zero)
%---------------------------------------------------------------------------
p  = p - OffSet;

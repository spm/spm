function f = spm_invFcdf(p,df)
% Inverse Cumulative Distribution Function (CDF) of F-distribution
% FORMAT f = spm_invFcdf(p,df)
% p   - Lower tail probability
% df  - degrees of freedom
%
%     - Multiple p may be supplied, with either constant df, or matrix of
%       pairs of df corresponding to the elements of p(:)'
%     - df must be a 2 x n (or n x 2) matrix of pairs of degrees of freedom.
%       (2x2 matrices are read as having df pairs in columns)
%__________________________________________________________________________
%
% spm_invFcdf implements the inverse Cumulative Distribution Function
% of F distributions.
%
% Returns the F-variate f such that Pr(F <= f) = p for F, a F random variable
% on df degrees of freedom.
%
% Works using spm_fzero to find the point where spm_Fcdf equals p.
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version control-%
% V1a	- 14/12/93 - 
% V1b	- 24/08/94 - streamlined code

%-Parameters
%---------------------------------------------------------------------------
Tol=[];

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end

if any(abs(p(:)-0.5)>0.5) error('p must be in [0,1]'), end
if any(df(:)<=0) error('df must be strictly positive'), end
if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end

%-re-orient df to 2 x n size, if n x 2 with n~=2.
if (size(df,1)~=2)&(size(df,2)==2) df=df'; end
if size(df,1)~=2 error('df must have 2 rows or 2 columns'), end

%-check sizes of arguments
if prod(size(p))==size(df,2)
elseif length(p)==1
	p=p*ones(1,size(df,2));
elseif size(df,2)==1
	df=meshgrid(df,1:prod(size(p)))';
else
	error('p and df not of compatible size'), end
end % if (size)

%-Computation
%---------------------------------------------------------------------------
f=zeros(size(p));
%-Avoid case where p=1, where invcdf is infinite
f(p==1)=+Inf*ones(sum(p==1),1);

InitGuess=10;
trace=0;
%-Avoid case where p=0, since Fcdf(f,df) is zero for all negative f
for k=find(abs(p-0.5)<0.5)
	f(k)=spm_fzero('spm_Fcdf',InitGuess,Tol,trace,df(:,k),p(k));
end % for

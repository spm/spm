function pdf = spm_Fpdf(f,df)
% Probability Density Function (PDF) of F distribution
% FORMAT pdf = spm_Fpdf(f,df)
% f      - F-variate
% df     - (vector of) degrees of freedom
%
%        - Multiple f may be supplied, with either constant df, or matrix of
%          pairs of df corresponding to the elements of f(:)'
%        - df must be a 2 x n (or n x 2) matrix of pairs of degrees of 
%          freedom. (2x2 matrices are read as having df pairs in columns)
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version control-%
% V1b	- 24/08/94 - streamlined code

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end

if any(df(:)<=0) error('df must be strictly positive'), end
if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end

%-re-orient df to 2 x n size, if n x 2 with n~=2.
if (size(df,1)~=2)&(size(df,2)==2) df=df'; end
if size(df,1)~=2 error('df must have 2 rows or 2 columns'), end

%-check sizes of arguments
if prod(size(f))==size(df,2)
elseif length(f)==1
	f=f*ones(1,size(df,2));
elseif size(df,2)==1
	df=meshgrid(df,1:prod(size(f)))';
else
	error('f and df not of compatible size'), end
end % if (size)

%-Store sizes-reshape arguments to rows
%---------------------------------------------------------------------------
[rdim,cdim]=size(f);
f=f(:)';
n=length(f);

%-Computation
%---------------------------------------------------------------------------
%-avoid cases f<=0 where pdf=0
pdf=zeros(n,1);
K=find(f>0);
if length(K)>0;
	df1K=df(1,K); df2K=df(2,K);
	scK=1./beta(df2K/2,df1K/2);
	pdf(K)=scK.*(df1K./df2K).^(df1K/2).*f(K).^(df1K/2-1).*...
		(1+(df1K./df2K).*f(K)).^(-(df1K+df2K)/2);
end % f>0 cases

pdf=reshape(pdf,rdim,cdim);

function p = spm_Gcdf(g,df,OffSet)
% Cumulative Distribution Function (CDF) of Gamma distribution
% FORMAT p = spm_Gcdf(g,df[,OffSet])
% g      - Gamma-variate
% df     - (vector of) degrees of freedom
% OffSet - (optional) offset, subtracted from resulting probability.
%          This enables use with fzero to compute inverse CDF.
%
%        - Multiple g may be supplied, with either constant df, or matrix of
%          pairs of df corresponding to the elements of g(:)'
%        - df must be a 2 x n (or n x 2) matrix of pairs of degrees of
%          freedom. (2x2 matrices are read as having df pairs in columns)
%__________________________________________________________________________
%
% spm_Gcdf implements the Cumulative Distribution of the Gamma-distribution.
%
% Returns the probability p, that a Gamma-variate on df degrees of freedom
% is less than g. p = Pr(G <= g) for G, a Gamma random variable on df degrees 
% of freedom.
%
% Uses the incomplete gamma function, as described in
% "Numerical Recipies in C"
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<3 OffSet=0; end
if nargin<2 error('insufficient arguments'), end

if any(abs(OffSet(:)-0.5)>0.5) error('OffSet must be in [0,1]'), end
if any(df(:)<=0) error('df must be strictly positive'), end

%-re-orient df to 2 x n size, if n x 2 with n~=2.
if (size(df,1)~=2)&(size(df,2)==2) df=df'; end
if size(df,1)~=2 error('df must have 2 rows or 2 columns'), end

%-check sizes of arguments
if prod(size(g))==size(df,2)
elseif length(g)==1
	g=g*ones(1,size(df,2));
elseif size(df,2)==1
	df=meshgrid(df,1:prod(size(g)))';
else
	error('g and df not of compatible size'), end
end % if (size)

%-Store sizes-reshape arguments to rows
%---------------------------------------------------------------------------
[p_rdim,p_cdim]=size(g);
g=g(:)';
n=length(g);

%-Computation
%---------------------------------------------------------------------------
%-avoid cases g<=0 where p=0
p=zeros(1,n);
K=find(g>0);
if length(K)>0;
	df1K=df(1,K); df2K=df(2,K);
	p(K)=gammainc(df2K.*g(K),df1K);
end % g>0 cases

p=reshape(p,p_rdim,p_cdim); end

%-Subtract (which defaults to zero)
%---------------------------------------------------------------------------
p=p-OffSet;

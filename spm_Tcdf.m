function p = spm_Tcdf(t,df,OffSet)
% Cumulative Distribution Function (CDF) of Students t distribution
% FORMAT p = spm_Tcdf(t,df[,OffSet])
% t      - t-variate
% df     - degrees of freedom
%        - t and df must be compatible for addition.
% OffSet - (optional) offset, subtracted from resulting probability.
%          This enables use with fzero to compute inverse CDF.
%__________________________________________________________________________
%
% spm_Tcdf implements the Cumulative Distribution of Students t-distributions.
%
% Returns the probability p, that a Students t-variate on df degrees of 
% freedom is less than t. p = Pr(T <= t) for T, a Students t random
% variable on df degrees of freedom.
%
% Uses the incomplete beta function, as described in "Numerical Recipies in C"
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version control-%
% V2a	- 06/08/93 - used numerical integration 
% V2b	- 14/12/93 - use incomplete Beta function
% V2c	- 24/08/94 - streamlined code

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<3 OffSet=0; end
if nargin<2 error('insufficient arguments'), end

if any(abs(OffSet(:)-0.5)>0.5) error('OffSet must be in [0,1]'), end
if any(df(:)<=0) error('df must be strictly positive'), end
if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end

%-if 1scalar & 1matrix argument, extend scalar one to matrix size
if all(size(t)==size(df))
elseif length(df)==1
	df=df*ones(size(t));
elseif length(t)==1
	t=t*ones(size(df));
else
	error('t and df not compatible for addition'), end
end % if (size)

p=zeros(size(t));

%-Computation
%---------------------------------------------------------------------------
%-avoid case t==0 to avoid log of zero warning, t=0 => p=0.5
K=find(t==0);
if length(K)>0
	p(K)=0.5 * ones(length(K),1);
end % t==0 cases

K=find(t~=0);
if length(K)>0;
	p(K)=0.5*betainc(df(K)./(df(K)+t(K).^2),df(K)/2,1/2);
	p(t>0)=1-p(t>0);
end % t~=0 cases

%-Subtract (which defaults to zero)
%---------------------------------------------------------------------------
p=p-OffSet;

function t = spm_invTcdf(p,df)
% Inverse Cumulative Distribution Function (CDF) of Students t distribution
% FORMAT t = spm_invTcdf(p,df)
% p   - Lower tail probability
% df  - degrees of freedom
%     - p and df must be compatible for addition.
%__________________________________________________________________________
%
% spm_invTcdf implements the inverse Cumulative Distribution of Students
% t-distributions.
%
% Returns the t-variate t such that Pr(T <= t) = p for T, a Students t
% random variable on df degrees of freedom.
%
% t=spm_invTcdf(1-p,df) therefore returns the t-variate corresponding to a
% (probability) threshold of p : Pr(T>t)=p.
%
% Works using spm_fzero to find the point where spm_Tcdf equals p.
%
% p &/or df can be matricies, in which case elements of p and df are
% paired. p & df must be compatible for addition.
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-Parameters
%---------------------------------------------------------------------------
Tol=[];

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end

if any(abs(p(:)-0.5)>0.5) error('p must be in [0,1]'), end
if any(df(:)<=0) error('df must be strictly positive'), end
% if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end

%-Check sizes of arguments
%-if 1scalar & 1matrix argument, extend scalar one to matrix size
if all(size(p)==size(df))
elseif length(df)==1
	df=df*ones(size(p));
elseif length(p)==1
	p=p*ones(size(df));
else
	error('p and df not compatible for addition')
end

%-Computation
%---------------------------------------------------------------------------
t = zeros(size(p));
%-Avoid cases where p=0,1, where invcdf is infinite
t(p==0)=-Inf*ones(sum(p==0),1);
t(p==1)=+Inf*ones(sum(p==1),1);

InitGuess=0;
trace=0;
for k=find(abs(p-0.5)<0.5)
	t(k)=spm_fzero('spm_Tcdf',InitGuess,Tol,trace,df(k),p(k));
end

function pdf = spm_Gpdf(g,df)
% Probability Density Function (PDF) of Gamma distribution
% FORMAT pdf = spm_Gpdf(g,df)
% g      - Gamma-variate
% df     - (vector of) degrees of freedom
% pdf    - PDF of Gamma-distribution with df degrees of freedom at points t
%
%        - Multiple g may be supplied, with either constant df, or matrix of
%          pairs of df corresponding to the elements of g(:)'
%        - df must be a 2 x n (or n x 2) matrix of pairs of degrees of 
%          freedom. (2x2 matrices are read as having df pairs in columns)
%__________________________________________________________________________
%
% spm_Gpdf implements the Probability Density Function of the Gamma
% distribution.
%
%---------------------------------------------------------------------------
% %W% Andrew Holmes %E%

%-version control-%
% V1a	- 13/12/93 - 
% V1b	- 24/08/94 - streamlined code

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end

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
[rdim,cdim]=size(g);
g=g(:)';
n=length(g);

%-Computation
%---------------------------------------------------------------------------
%-avoid cases g<=0 where pdf=0
pdf=zeros(1,n);
K=find(g>0);
if length(K)>0;
	df1K=df(1,K); df2K=df(2,K);
	pdf(K)=(df2K).^(df1K).*g(K).^(df1K-1).*exp(-df2K.*g(K))...
					./gamma(df1K);
end % g~=0 cases

pdf=reshape(pdf,rdim,cdim);

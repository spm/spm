function p = spm_Xcdf(x,df,OffSet)
% Cumulative Distribution Function (CDF) of Chi-squared distribution
% FORMAT p = spm_Xcdf(x,df,OffSet)
% x      - Chi-squared variate
% df     - degrees of freedom
%        - x and df must be compatible for addition.
% OffSet - (optional) offset, subtracted from resulting probability.
%          This enables use with fzero to compute inverse CDF.
%__________________________________________________________________________
%
%
% spm_Xcdf implements the Cumulative Distribution of Chi-squared
% distributions.
%
% Returns the probability p, that a Students Chi-squared variate on df
% degrees of freedom is less than x. p = Pr(X <= x) for X a Chi-squared
% random variable on df degrees of freedom.
%
% A Chi-Squared distribution on df degrees of freedom is a Gamma
% distribution with [df/2,1/2] degrees of freedom. This identity
% is used to compute the PDF with spm_Gcdf
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<3 OffSet=0; end
if nargin<2 error('insufficient arguments'), end

if any(abs(OffSet(:)-0.5)>0.5) error('OffSet must be in [0,1]'), end
if any(df(:)<=0) error('df must be strictly positive'), end
% if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end

%-Computation
%---------------------------------------------------------------------------
Gdf=[df(:)'/2; ones(1,length(df(:)))/2];

p=spm_Gcdf(x,Gdf);

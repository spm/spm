function x = spm_invXcdf(p,df)
% Inverse Cumulative Distribution Function (CDF) of Chi-squared distribution
% FORMAT x = spm_invXcdf(p,df)
% p   - Lower tail probability
% df  - degrees of freedom
%     - p and df must be compatible for addition.
%__________________________________________________________________________
%
% spm_invXcdf implements the inverse Cumulative Distribution of the
% Chi-squared distributions.
%
% Returns the Chi-squared variate x such that Pr(X <= x) = p for X a 
% Chi-squared random variable on df degrees of freedom.
%
% x = spm_invTcdf(1 - p,df) therefore returns the Chi-squared variate 
% corresponding to a (probability) threshold of p : Pr(T > t) = p.
%
% Works using spm_fzero to find the point where spm_Xcdf equals p.
%
% A Chi-Squared distribution on df degrees of freedom is a Gamma
% distribution with [df/2,1/2] degrees of freedom. This identity
% is used to compute the PDF with spm_Gcdf
%
% p &/or df can be matrices, in which case elements of p and df are
% paired. p & df must be compatible for addition.
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version control-%
% V1a	- 14/12/93 - Andrew Holmes

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end

if any(abs(p(:)-0.5)>0.5) error('p must be in [0,1]'), end
if any(df(:)<=0) error('df out of range'), end

%-Computation
%---------------------------------------------------------------------------
Gdf=[df(:)'/2; ones(1,length(df(:)))/2];

x=spm_invGcdf(p,Gdf);

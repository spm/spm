function pdf = spm_Xpdf(x,df)
% Probability Density Function (PDF) of Chi-Squared distribution
% FORMAT pdf = spm_Xpdf(x,df)
% x      - Chi-squared variate
% df     - degrees of freedom
%        - x and df must be compatible for addition.
% pdf    - PDF of Chi-squared distribution with df degrees of freedom,
%          at points x
%__________________________________________________________________________
%
% spm_Xpdf implements the Probability Density Function of the 
% Chi-squared distributions.
%
% A Chi-Squared distribution on df degrees of freedom is a Gamma
% distribution with [df/2,1/2] degrees of freedom. This identity
% is used to compute the PDF with spm_Gpdf
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version control-%
% V1a	- 13/12/93 - Andrew Holmes

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end

if any(df(:)<=0) error('df must be strictly positive'), end
if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end

%-Computation
%---------------------------------------------------------------------------
Gdf=[df(:)'/2; ones(1,length(df(:)))/2];

pdf=spm_Gpdf(x,Gdf);

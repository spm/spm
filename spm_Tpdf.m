function pdf = spm_Tpdf(t,df)
% Probability Density Function (PDF) of Students t distribution
% FORMAT pdf = spm_Tpdf(t,df)
% t      - t-variate
% df     - degrees of freedom
%        - t and df must be compatible for addition.
% pdf    - PDF of t-distribution with df degrees of freedom at points t
%__________________________________________________________________________
%
% spm_Tpdf implements the Probability Density Function of Students 
% t-distributions.
%
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-version control-%
% V2a	- 06/08/93 - used two gamma functions
% V2b	- 13/12/93 - changed to use one Beta function call

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end

if any(df(:)<=0) error('df must be strictly positive'), end
if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end

%-Computation
%---------------------------------------------------------------------------
pdf=((1+t.^2./df).^(-(df+1)/2)) ./ (sqrt(df).*beta(df/2,1/2));

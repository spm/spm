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

%-Argument range and size checks
%---------------------------------------------------------------------------
if nargin<2 error('insufficient arguments'), end
% if any(floor(df(:))~=ceil(df(:))) error('df must be integer'), end
if any(df(:)<=0) error('df must be strictly positive'), end

%-Computation
%---------------------------------------------------------------------------
pdf=((1+t.^2./df).^(-(df+1)/2)) ./ (sqrt(df).*beta(df/2,1/2));

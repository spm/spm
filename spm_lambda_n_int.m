function y = spm_lambda_n_int(t,n)
% Integrand of Lambda_n, relating smoothness of Gaussianised and raw t-fields
% FORMAT y = spm_lambda_n_int(t,n)
%
% t - t variable
% n - df + 1 : Degrees of freedon of the t-field, plus one
%_______________________________________________________________________
%
% Integrand of \lambda_n, the ratio of the Variance-Covariance matrix
% of partial derivatives of a Gaussianised t-field to the components of
% the t-field.
% 
% \lambda_n is the integral of spm_lambda_n_int over all real values.
%
% The formula is that on p917 of Worsley et al. (1992) Worsley et al.,
% (1992) "A Three Dimensional Statistical Analysis for CBF Activation
% Studies in Human Brain", Journal of Cerebral Blood and Metabolism
% 12:900-918
%
%_______________________________________________________________________
% %W% Stefan Kiebel, Andrew Holmes %E%

% Check arguments
%-----------------------------------------------------------------------
if (nargin < 2), y = []; return, end

% Computation
%=======================================================================
tmp1 = (t.^2 + n - 1).^2/((n-1)*(n-2));
tmp2 = (spm_Tpdf(t, n - 1)).^3;
tmp3 = (spm_Npdf(spm_invNcdf(1 - spm_Tcdf(t, n - 1)))).^2;
y = tmp1.*tmp2./tmp3;

return

% Notes
%=======================================================================
% tmp1 = (t.^2 + n - 1).^2;
% tmp2 = (n-1)*(n-2);
% tmp3 = spm_Tpdf(t, n-1);
% tmp3 = (spm_Npdf(spm_invNcdf(spm_Tcdf(1 - t, n - 1)))).^2;
% value = tmp1.*tmp3./tmp2;

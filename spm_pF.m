function CPF = spm_pF(S,W,df,F)
% Returns the corrected p value for a SPM{F}
% FORMAT CPF = spm_pF(S,W,df,F)
% F   - value of the F-field
% S   - volume analysed
% W   - smoothness of the component fields
% df  - [df(interest) df(residuals)]
%___________________________________________________________________________
%
% Returns the corrected p value for the F statistic
% Reference : Worsley 1994, Adv Appl Prob 26, 13-42.
%
%---------------------------------------------------------------------------
% %W% JBP %E%

% check on df
%---------------------------------------------------------------------------
n      = df(1);
m      = df(2);
if (size(W,2) == 3 & (m < 4)) | (size(W,2) == 2 & (m < 3))
	error('df too small'); end

if size(W,1) > 1
	W   = mean(W); end

% parameters
%---------------------------------------------------------------------------
SqrtDetLmda = prod(1./(W*sqrt(2)));
Lm2 	    = gamma(m/2);
Ln2 	    = gamma(n/2);
f	    = n*F/m;

% compute expection {2D}
%---------------------------------------------------------------------------
if size(W,2) == 2
	CPF	= S*SqrtDetLmda*gamma((m + n-2)/2)*...
	(f^(n/2 -1)) * ((1 + f)^(-(m + n - 2)/2)) * ((m - 1)*f - n + 1)/...
	(2*pi*Lm2*Ln2);
end

% compute expection {3D}
%---------------------------------------------------------------------------
if size(W,2) == 3
	CPF	= S*SqrtDetLmda*gamma((m + n - 3)/2)*...
	f^(n/2 - 3/2) * ((1+f)^(-(m+n-2)/2) ) * ...
	((m - 1)*(m - 2)*(f^2) - (2*m*n - m - n - 1)*f + (n - 1)*(n - 2))/...
	((2*pi)^(3/2)*sqrt(2)*Lm2*Ln2);

end

% p value assuming Poisson behaviour
%---------------------------------------------------------------------------
CPF = 1 - exp(-CPF); 



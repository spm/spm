function CPF = spm_pF(S, Wresid, Fdf, Fval)

% function CPF = spm_pF(S, Wresid, Fdf, Fval)
%
%
% Return the corrected p value for the F statistic
% Reference : Worsley 1994, Adv Appl Prob 26, 13-42.
%
% S : volume analysed
% Wresid : smoothness of the component fields
% Fdf : [df(of interest) df(residuals)]
% Fval : value of the F-field
% %W% JBP %E%


%---------------------------------------------------------------------------
n      = Fdf(1); m      = Fdf(2);
if (size(Wresid,2) == 3 & m<4) | (size(Wresid,2) == 2 & m<3)
	disp(['ERROR :  df too small']); return; end

%---------------------------------------------------------------------------
SqrtDetLmda 	= sqrt(det(diag(mean(((Wresid.^2)*2).^(-1)))));
Lm2 		= gamma(m/2);
Ln2 		= gamma(n/2);
f		= n*Fval/m;

if size(Wresid,2) == 2
   CPF	= S*SqrtDetLmda*gamma((m+n-2)/2)*...
		(f^(n/2 -1)) * ((1+f)^(-(m+n-2)/2)) * ((m-1)*f-n+1) /...
		( 2*pi*Lm2*Ln2 );
end

if size(Wresid,2) == 3
   CPF	= S*SqrtDetLmda*gamma((m+n-3)/2)*...
		f^(n/2 - 3/2) * ((1+f)^(-(m+n-2)/2) ) * ...
		( (m-1)*(m-2)*(f^2) - (2*m*n - m - n - 1)*f +(n-1)*(n-2) ) /...
		( (2*pi)^(3/2) * sqrt(2) * Lm2 * Ln2 );

end
CPF = 1 - exp(-CPF); % Suppose Poisson behaviour



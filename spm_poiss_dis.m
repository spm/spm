function [y] = spm_poiss_dis(x,q)
% returns the t distribution
% FORMAT [y] = spm_poiss_dis(x,q)
% x   - range in real positive integers
% q   - parameter
%___________________________________________________________________________
%
% spm_poiss_dis return the Poisson distribution with parameter q
% over the range x
%
%__________________________________________________________________________
% %W% %E%

fprintf('%s\nWARNING: spm_poiss_dis is grandfathered: Use spm_Ppdf instead...\n',7)

y = spm_Ppdf(x,q);

return




%---------------------------------------------------------------------------
X = x;
x = 1:max(x);

if q > 64
	y = [exp(-q) 1/sqrt(2*pi*q)*exp(-(x - q).^2/(2*q))];
else
	y = [exp(-q) exp(-q)*q.^x./cumprod(x)];
end

y = y(X + 1);

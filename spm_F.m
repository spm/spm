function [u] = spm_F(a,s,Fdf,S)
% critical height threshold at a specified significance level
% FORMAT [u] = spm_F(a,s,Fdf,S)
% u   - critical height spm_pF(S,s,Fdf,u) = a
%
% a   - level of significance - alpha (eg 0.05)
% s   - smoothness - length(s) = D - dimension
% Fdf - df of F distribution
% S   - Lebesgue measure of S
%___________________________________________________________________________
% spm_F returns the critical threshold at a specified significance
% volume S of a D-dimensional SPM{F} with component smoothness s.
%
%  Reference : Worsley 1994, Adv Appl Prob 26, 13-42.
%___________________________________________________________________________
% %W% Karl Friston %E%

% Gauss-Newton search
%---------------------------------------------------------------------------
RESEL = S/prod(sqrt(8*log(2))*s);
u     = spm_invFcdf((1 - a/RESEL),Fdf);
du    = 1e-4;

% refined estimate using 1 - exp(-E{m})
%---------------------------------------------------------------------------
d     = 1;
while abs(d) > 1e-3
	q     = spm_pF(S,s,Fdf,u);
	dqdu  = (spm_pF(S,s,Fdf,u + du) - q)/du;
	d     = (q - a)/dqdu;
	u     = u - d;
end

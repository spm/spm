function [u] = spm_z(a,s,S)
% critical height threshold at a specified significance level
% FORMAT [u] = spm_z(a,s,S)
% u   - critical height P(u,k,c) = a
%
% a   - level of significance - alpha (eg 0.05)
% s   - smoothness - length(s) = D - dimension
% S   - Lebesgue measure of S
%___________________________________________________________________________
% spm_u returns the critical threshold at a specified significance
% volume S of a D-dimensional Gaussian process of isotropic smoothness s.
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% Ref: Friston et al (1994) Assessing the significance of focal
% activations using their spatial extent. Human Brain Mapping 
% 1:210-220

%___________________________________________________________________________
% %W% Karl Friston %E%

% Gauss-Newton search
%---------------------------------------------------------------------------
D     = length(s);
s     = prod(s).^(1/D);
d     = 1;
u     = 4;
du    = 1e-4;

% find approximate value using E{m}
%---------------------------------------------------------------------------
while abs(d) > 1e-3
	p     = S*(2*pi)^(-(D + 1)/2)*(2*s^2)^(-D/2)*u^(D - 1)*exp(-(u^2)/2);
	u     = u + du;
	q     = S*(2*pi)^(-(D + 1)/2)*(2*s^2)^(-D/2)*u^(D - 1)*exp(-(u^2)/2);

	dEdu  = (q - p)/du;
	d     = (q - a)/dEdu;
	u     = u - d;
end

% refined estimate using 1 - exp(-E{m})
%---------------------------------------------------------------------------
d     = 1;
while d > 1e-3
	E     = S*(2*pi)^(-(D + 1)/2)*(2*s^2)^(-D/2)*u^(D - 1)*exp(-(u^2)/2);
	p     = 1 - exp(-E);

	u     = u + du;
	H     = S*(2*pi)^(-(D + 1)/2)*(2*s^2)^(-D/2)*u^(D - 1)*exp(-(u^2)/2);
	q     = 1 - exp(-H);

	dEdu  = (q - p)/du;
	d     = (q - a)/dEdu;
	u     = u - d;
end

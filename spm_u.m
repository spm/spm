function [K] = spm_u(a,s,u,D,S)
% critical resgion size at a specified significance level
% FORMAT [u] = spm_u(a,s,u,D,S)
% u   - critical size P(nmax >= u) = alpha
% a   - level of significance - alpha (eg 0.05)
% s   - smoothness
% u   - threshold
% D   - dimension
% S   - Lebesgue measure of S
%___________________________________________________________________________
% spm_u returns the  critical region size at a specified significance
% and threshold in volume S of a D-dimensional Gaussian process
% of isotropic smoothness s, thresholded at u.
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% Ref: Friston et al (1993) Comparing functional images: Assessing
% the spatial extent of activation foci
%
%__________________________________________________________________________
% %W% %E%


%---------------------------------------------------------------------------
EN       = S*(1 - spm_Ncdf(u));
Em       = S*(2*pi)^(-(D + 1)/2)*(2*s^2)^(-D/2)*u^(D - 1)*exp(-(u^2)/2);
En       = EN/Em;
b        = (gamma(D/2 + 1)/En)^(2/D);
K        = (-log(-log(1 - a)/Em)/b)^(D/2);

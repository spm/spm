function [P] = spm_Pn(s,u,S)
% probability based on maximum P(Zmax >= u)
% FORMAT [P] = spm_Pn(s,u,S);
% s   - smoothness: length(s) =  D - dimension
% u   - threshold
% S   - Lebesgue measure of S {volume in voxels}
% p   - P(Zmax >= u)
%___________________________________________________________________________
%
% Returns the probability of one or more regions with a height
% greater than u in volume S of a D - dimensional Gaussian process
% of non-isotropic smoothness s
%
% The value returned uses the first moment [expectation <m>] of the p.d.f
% of the number of maxima m to estimate P(m > 0) = 1 - exp(-<m>) assuming
% a Poisson form for P(m = M)
%
% see also spm_Pn.m
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
%
%___________________________________________________________________________
% %W% %E%

%---------------------------------------------------------------------------
D = length(s);
P = S*(2*pi)^(-(D + 1)/2)*prod(2*s.^2)^(-1/2)*u.^(D - 1).*exp(-(u.^2)/2);
P = 1 - exp(-P);


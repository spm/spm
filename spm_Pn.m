function [P] = spm_Pn(k,s,u,S)
% probility based on number of voxels P(nmax >= k)
% FORMAT [P] = spm_Pn(k,s,u,S);
% k   - number of voxels in the region
% s   - smoothness: length(s) =  D = dimension
% u   - threshold
% S   - Lebesgue measure of S {volume in voxels}
% p   - P(nmax >= k)
%___________________________________________________________________________
%
% spm_Pn returns the probability of one or more regions with more than
% k voxels in volume S of a D - dimensional Gaussian process of
% non-isotropic smoothness s, thresholded at u.
%
% see also spm_Pu.m
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% Ref: Friston et al (1994) Assessing the significance of focal
% activations with their spatial extent. HBM 1:210-220
%
%__________________________________________________________________________
% %W% %E%

% EN = expected total number of voxels {N}
% Em = expected total number of maxima {m}
% Em = expected total number of voxels per region (or maxima) {n}
% b  = parameter for the pdf P(n = X) {See Friston et al}

%---------------------------------------------------------------------------
D    = length(s);
EN   = S*(1-spm_Ncdf(u));
Em   = S*(2*pi)^(-(D + 1)/2)*prod(2*s.^2)^(-1/2)*u^(D - 1)*exp(-(u^2)/2);
En   = EN/Em;
b    = (gamma(D/2 + 1)/En)^(2/D);
P    = 1 - exp(-Em*exp(-b*k.^(2/D)));

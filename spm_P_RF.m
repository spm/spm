function [P,EN,Em,En,Pk] = spm_P(c,s,u,k,S)
% probability based on the number of clusters P(C > c)
% FORMAT [P,EN,Em,En,Pk] = spm_P(c,s,u,k,S);
%
% c   - number of voxels in the region
% s   - smoothness: length(s) = D {dimension}
% u   - height threshold
% k   - extent threshold
% S   - Lebesgue measure of S {volume in voxels}
%
% P   - P(u,k,c)
% EN  - expected total number of voxels {N}
% Em  - expected total number of maxima {m}
% En  - expected total number of voxels per region (or maxima) {n}
% b   - parameter for the pdf P(n = x)
% Pk  - P(n > k)
%___________________________________________________________________________
%
% spm_P returns the probability of c or more regions with more than
% k voxels in volume S of a D - dimensional Gaussian process of
% non-isotropic smoothness s, thresholded at u.
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% Ref: Friston et al (1994) Assessing the significance of focal
% activations with their spatial extent. HBM 1:210-220
%
%__________________________________________________________________________
% %W% Karl Friston %E%


%---------------------------------------------------------------------------
D    = length(s);
u    = min([u 24]);
if u > 6
	EN = S*exp(-(u.^2)/2)./(u*sqrt(2*pi));
else
	EN = S*(1 - spm_Ncdf(u));
end
Em   = S*(2*pi)^(-(D + 1)/2)*prod(2*s.^2)^(-1/2)*u^(D - 1)*exp(-(u^2)/2);
En   = EN/Em;
beta = (gamma(D/2 + 1)/En)^(2/D);
Pk   = exp(-beta*(k^(2/D)));
P    = 1 - sum(spm_Ppdf([0:(c - 1)],Em*Pk));

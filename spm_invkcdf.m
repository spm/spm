function k = spm_invkcdf(p,u,s)
% Inverse CDF for k, the number of voxels per cluster
% FORMAT k = spm_invkcdf(p,u,s);
%
% x   - lower tail probability
% s   - smoothness: length(s) = D {dimension}
% u   - height threshold
% 
% k   - size variate
%__________________________________________________________________________
%
% spm_kcdf implements the inverse of Cumulative Distribution Function (CDF)
% for the number of voxels in a given cluster, conditional upon that cluster 
% existing.
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% Ref: Friston et al (1994) Assessing the significance of focal
% activations with their spatial extent. HBM 1:210-220
%
%__________________________________________________________________________
% %W% Karl Friston %E%

% EN   - expected total number of voxels {N}
% Em   - expected total number of maxima {m}
% En   - expected total number of voxels per region (or maxima) {n}
% beta - parameter for the pdf P(n = x)

% determine CDF
%---------------------------------------------------------------------------
D       = length(s);
EN      = (1 - spm_Ncdf(u));
Em      = (2*pi)^(-(D + 1)/2)*prod(2*s.^2)^(-1/2)*u^(D - 1)*exp(-(u^2)/2);
En      = EN/Em;
beta    = (gamma(D/2 + 1)/En)^(2/D);
k       = (-log(1 - p)/beta).^(D/2);

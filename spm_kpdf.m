function pdf = spm_kcdf(x,u,s)
% PDF for k, the number of voxels per cluster
% FORMAT pdf = spm_kcdf(x,u,s);
%
% x   - ordinate - number of voxels in the region
% s   - smoothness: length(s) = D {dimension}
% u   - height threshold
% 
% cdf - P(k < x);
%__________________________________________________________________________
%
% spm_kpdf implements the Probability Density Function (PDF) for the
% number of voxels in a given cluster, conditional upon that cluster existing.
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

% determine PDF
%---------------------------------------------------------------------------
D       = length(s);
EN      = (1 - spm_Ncdf(u));
Em      = (2*pi)^(-(D + 1)/2)*prod(2*s.^2)^(-1/2)*u^(D - 1)*exp(-(u^2)/2);
En      = EN/Em;
beta    = (gamma(D/2 + 1)/En)^(2/D);
pdf     = (2*beta/D)*(x.^(2/D - 1)).*exp(-beta*(x.^(2/D)));

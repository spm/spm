function [S] = spm_fit_fmin(P,x,y)
% minimization function called by spm_fit_Gaussians
% FORMAT [S] = spm_fit_fmin(p,x,y);
% p[1] : mean 1
% p[2] : sigma 1
% p[3] : mean 2
% p[4] : sigma 2
% x    - bins
% y    - frequency
% S    - sum of squares
%___________________________________________________________________________
%
% spm_fit_fmin is the objective function minimized by spm_fit_Gaussians.
% The two Gaussians are contrained to have equal integrals.  This generally
% gives better convergence and a more acceptable solution when the fit is
% poor.
%
%__________________________________________________________________________
% %W% %E%

% ratio of integrals = r
%---------------------------------------------------------------------------
r       = 0.5;
g1      = exp(-(x - P(1)).^2/(2*P(2)^2));
g1      = r*g1/sum(g1);
g2      = exp(-(x - P(3)).^2/(2*P(4)^2));
g2      = (1 - r)*g2/sum(g2);

% sum of squares
%---------------------------------------------------------------------------
S       = sum((g1 + g2 - y/sum(y)).^2);



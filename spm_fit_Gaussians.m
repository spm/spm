function [p] = spm_fit_Gaussians(x,y)
% fits two Gaussians to a histogram
% FORMAT [p] = spm_fit_Gaussians(x,y);
% x     - Bins
% y     - Frequency
% p[1] : mean 1
% p[2] : sigma 1
% p[3] : mean 2
% p[4] : sigma 2
%___________________________________________________________________________
%
% spm_fit_Gaussians fits two Gaussians to a histogram using nonlinear
% minimization.  The two Gaussians are constrained to have integrals
% that bear some fixed relation (specified in the objective function
% - spm_fit_fmin.m).
%
% This routine is used primarily to estimate segmentation functions
% in application to MRI data.
%
%__________________________________________________________________________
% %W% %E%

% initial parameters
%---------------------------------------------------------------------------
maxx    = max(x);
y       = y.*(y > max(y)/5);
x0(1)   = 0.5*maxx;
x0(2)   = 0.1*maxx;
x0(3)   = 0.6*maxx;
x0(4)   = 0.1*maxx;

% nonlinear minimization in parameter space
%---------------------------------------------------------------------------
p       = fmins('spm_fit_fmin',x0,[0, 0.0001],[],x(:),y(:));


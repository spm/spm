function [DCM] = spm_dcm_fmri_check(P)
% post-hoc diagnostics for DCM (bilinear or nonlinear) of fMRI data
% FORMAT [DCM] = spm_dcm_fmri_check(DCM)
%   DCM - DCM structure or its filename
%
% This routine provides some diagnostics to ensure model inversion has
% converged. It plots the predicted and observed responses over all regions
% and provides the coefficient of determination - or percent variance
% explained. This should normally be above 10%. An abnormally low
% coefficient of determination is highlighted in red. Quantitatively, one
% would normally expect to see one or more extrinsic (between source)
% connections with the strength of 1/8 Hz or greater. If all the extrinsic
% posterior expectations are below this value, then this suggests a failure
% of convergence or that the data are very noisy (possibly due to using
% very small regions of interest to summarise regional responses). Finally,
% the posterior correlations among all parameters are shown in terms of a
% correlation matrix. The number of effective parameters estimated is
% reported in terms of the (KL) divergence between the posterior and
% prior densities over parameters. This is divided by the log of the
% number of observations, by appealing to the Bayesian information
% criterion. The divergence corresponds to complexity or Bayesian
% surprise. Normally, one would expect the posterior and prior to diverge
% in a non-trivial fashion.
%
% Posterior densities are shown as bars with 90% confidence intervals in
% pink. An informed model inversion would normally provide posterior
% densities with confidence intervals that are, for some connections,
% displaced from prior expectations (at or around zero).
%
% This routine is compatible with DCM8, DCM10 and DCM12 files.
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_check.m 5615 2013-08-15 14:37:24Z spm $


%-Load DCM structure
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, DCM = []; return; end
end
if isstruct(P)
    DCM = P;
else
    load(P)
end

% Assemble diagnostics
%==========================================================================

% coefficient of determination (percent variance explained)
%--------------------------------------------------------------------------
PSS   = sum(sum(DCM.y.^2));
RSS   = sum(sum(DCM.R.^2));
D(1)  = 100*PSS/(PSS + RSS);

% largest absolute posterior expectation (extrinsic connections)
%--------------------------------------------------------------------------
try
    A = DCM.Ep.A;
catch
    A = DCM.A;
end

if DCM.options.two_state
    A = exp(A);
end

D(2)  = max(max(abs(A - diag(diag(A)))));

% complexity and effective number of parameters estimated
%--------------------------------------------------------------------------
qE    = spm_vec(DCM.Ep);
pE    = spm_vec(DCM.M.pE);
qC    = DCM.Cp;
pC    = DCM.M.pC;
k     = rank(full(pC));
pC    = pinv(pC);

D(3)  = trace(pC*qC) + (pE - qE)'*pC*(pE - qE) - spm_logdet(qC*pC) - k;
D(3)  = D(3)/log(DCM.v);


% Plot summary of inversion
%==========================================================================
spm_figure('GetWin','DCM diagnostics'); clf


% plot predicted and observed regional responses
%--------------------------------------------------------------------------
subplot(2,1,1);
t   = (1:DCM.v)*DCM.Y.dt;
D   = full(D);

plot(t,DCM.y,t,DCM.y + DCM.R,':');
str = sprintf('variance explained %0.0f%%', D(1));
str = {'responses and predictions',str};
if D(1) > 10
    title(str,'FontSize',16);
else
    title(str,'FontSize',16,'Color','r');
end
xlabel('time {seconds}');


% posterior densities over A parameters
%--------------------------------------------------------------------------
try
    i = spm_fieldindices(DCM.Ep,'A');
catch
    i = 1 + (1:DCM.n^2);
end
qE  = spm_vec(DCM.Ep);
qC  = DCM.Cp;

if DCM.options.two_state
    qE = exp(qE);
end

subplot(2,2,3)
spm_plot_ci(qE(i),qC(i,i)), hold on
str = sprintf('largest connection strength %0.2f', D(2));
str = {'intrinsic and extrinsic connections',str};
if D(2) > 1/8
    title(str,'FontSize',16);
else
    title(str,'FontSize',16,'Color','r');
end
xlabel('parameters');
axis square


% posterior correlations among all parameters
%--------------------------------------------------------------------------
subplot(2,2,4)
imagesc(spm_cov2corr(DCM.Cp))
title('posterior correlations','FontSize',16)
str = sprintf('estimable parameters %0.0f', D(3));
str = {'posterior correlations',str};
if D(3) > 1
    title(str,'FontSize',16);
else
    title(str,'FontSize',16,'Color','r');
end
axis square

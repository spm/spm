function P = spm_eeg_inv_vbecd_mom(P)
% Model inversion routine for ECDs using variational Bayesian approach
% FORMAT P = spm_eeg_inv_vbecd_mom(P)
%
% Input:
% structure P with fields:
%  forward      - structure containing the forward model, i.e. the "vol"
%                 and "sens" structure in a FT compatible format
%  bad          - list of bad channels, not to use.
%  y            - data vector
%
%  Niter        - maximum number of iterations
%  priors       - priors on parameters,  as filled in (and
%                 described) in spm_eeg_inv_vbecd_gui.m.
%
% Output:
% same structure with extra fields
%  init         - initial valuse used for mu_w/s
%  dF           - successive (relative) improvement of F
%  post         - posterior value of estimated parameters and ther variance
%  Fi           - successive values of F
%  F            - Free energy final value.
%
% Reference:
% Kiebel et al., Variational Bayesian inversion of the equivalent current
% dipole model in EEG/MEG., NeuroImage, 39:728-741, 2008
% (Although this algorithm uses a function for general Bayesian inversion
% of a non-linear model - see spm_nlsi_gn)
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_inv_vbecd_mom.m 7764 2020-01-02 15:34:03Z spm $


% unpack model, priors, data
%--------------------------------------------------------------------------

priormom = P.priors.mom;
Nd    = length(priormom);   % number sources= number of moments to estimate
priormomvar = P.priors.momvar; % diagonal on covariance matrix for moment for each source

dippos_ctf = P.dippos_ctf;  % dipole positions, fixed
dipor_ctf  = P.dipor_ctf;   % dipole orientations, fixed.


y    = P.y;
Y.y  = y;

U.u  = 1;

% set random moment vector, scaled by prior variances
%--------------------------------------------------------------------------
startmoments = priormom + randn(Nd,1) .* sqrt(priormomvar);

% get lead fields
%--------------------------------------------------------------------------
M.pE = priormom;           % prior parameter estimate
M.pC = diag(priormomvar);  % prior covariance estimate

M.hE = P.priors.hE;
M.hC = P.priors.hC;

M.IS = 'spm_eeg_wrap_momfit_vbecd';

M.Setup = P;         % pass volume conductor and sensor locations on

M.pos = dippos_ctf;  % source positions
M.ori = dipor_ctf;   % source orientations


[starty] = spm_eeg_wrap_momfit_vbecd(startmoments,M,U);
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);
P.Ep = Ep;
P.Cp = Cp;
P.Eh = Eh;
P.F  = F;
[P.ypost,outsideflag,leads] = spm_eeg_wrap_momfit_vbecd(P.Ep,M,U);

P.post_mom = Ep;
P.post_momvar = Cp;

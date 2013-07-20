function DCM = spm_dcm_fmri_csd(DCM)
% Estimates parameters of a DCM (bilinear or nonlinear) for fMRI data
% FORMAT DCM = spm_dcm_fmri_csd(DCM)
%   DCM - DCM structure
%
% Expects
%--------------------------------------------------------------------------
% DCM.a                              % switch on endogenous connections
% DCM.b                              % switch on bilinear modulations
% DCM.c                              % switch on exogenous connections
% DCM.d                              % switch on nonlinear modulations
% DCM.U                              % exogenous inputs
% DCM.Y.y                            % responses
% DCM.Y.X0                           % confounds
% DCM.Y.Q                            % array of precision components
% DCM.n                              % number of regions
% DCM.v                              % number of scans
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_fmri_csd.m 5587 2013-07-20 15:37:17Z karl $
 

% check options
%==========================================================================
try, DCM.options.two_state;  catch, DCM.options.two_state  = 0;     end
try, DCM.options.stochastic; catch, DCM.options.stochastic = 0;     end
try, DCM.options.nonlinear;  catch, DCM.options.nonlinear  = 0;     end
try, DCM.options.centre;     catch, DCM.options.centre     = 0;     end
try, DCM.options.nmax;       catch, DCM.options.nmax       = 8;     end
try, DCM.options.hidden;     catch, DCM.options.hidden     = [];    end

try, M.nograph = DCM.options.nograph; catch, M.nograph = spm('CmdLine'); end
try, M.P       = DCM.options.P ;end

name = sprintf('DCM_%s',date);
DCM.options.analysis  = 'CSD';
 
% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                      catch, DCM.name = name;         end
try, DCM.options.Nmodes;            catch, DCM.options.Nmodes = 8;  end
 
% detrend outputs (and inputs)
%--------------------------------------------------------------------------
DCM.Y.y = spm_detrend(DCM.Y.y);
if DCM.options.centre
    DCM.U.u = spm_detrend(DCM.U.u);
end

% check scaling of Y (enforcing a maximum change of 4%
%--------------------------------------------------------------------------
scale       = max(max((DCM.Y.y))) - min(min((DCM.Y.y)));
scale       = 4/max(scale,4);
DCM.Y.y     = DCM.Y.y*scale;
DCM.Y.scale = scale;

% check confounds (add constant if necessary)
%--------------------------------------------------------------------------
if ~size(DCM.Y.X0,2), DCM.Y.X0 = ones(DCM.v,1); end


% check for endogenous DCMs, with no exogenous driving effects
%--------------------------------------------------------------------------
n    = DCM.n;
if isempty(DCM.c) || isempty(DCM.U.u)
    DCM.c      = zeros(n,1);
    DCM.b      = zeros(n,n,1);
    DCM.U.u    = zeros(v,1);
    DCM.U.name = {'null'};
end
if ~any(spm_vec(DCM.U.u)) || ~any(spm_vec(DCM.c))
    DCM.options.stochastic = 1;
end


% priors (and initial states)
%==========================================================================
[pE,pC,x] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);

% add prior on spectral density of fluctuations (amplitude and exponent)
%--------------------------------------------------------------------------
pE.a = sparse(2,n); pV.a = sparse(2,n) + 1/128; % neuronal fluctuations
pE.c = sparse(2,n); pV.c = sparse(2,n) + 1/128; % channel noise specific
pC   = spm_cat(spm_diag({pC,spm_diag(spm_vec(pV))}));




% create DCM
%--------------------------------------------------------------------------
DCM.M.IS = 'spm_csd_fmri_mtf';
DCM.M.FS = 'spm_fs_fmri_csd';
DCM.M.g  = 'spm_gx_fmri';
DCM.M.f  = 'spm_fx_fmri';
DCM.M.x  = x;
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.hE = 6;
DCM.M.hC = 1/128;
DCM.M.n  = length(spm_vec(x));
DCM.M.m  = size(DCM.U.u,2);
DCM.M.l  = n;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM.M.u   = sparse(n,1);

% get data-features
%==========================================================================
DCM       = spm_dcm_fmri_csd_data(DCM);
DCM.M.Hz  = DCM.Y.Hz;

% scale input (to a variance of one)
%--------------------------------------------------------------------------
ccf         = spm_csd2ccf(DCM.U.csd,DCM.Y.Hz);
DCM.U.scale = max(spm_vec(ccf));
DCM.U.csd   = spm_unvec(spm_vec(DCM.U.csd)/DCM.U.scale,(DCM.U.csd));


% complete model specification and invert
%==========================================================================

 
% precision of noise: AR(1/2)
%--------------------------------------------------------------------------
y        = spm_fs_fmri_csd(DCM.Y.csd,DCM.M);
n        = size(y,1);
m        = size(y,2)*size(y,3);
q        = spm_Q(1/2,n,1);
Q        = kron(speye(m,m),q);
DCM.Y.Q  = Q;
DCM.Y.X0 = sparse(size(Q,1),0);

spm_csd_fmri_mtf(pE,DCM.M,DCM.U);

% Variational Laplace: model inversion
%==========================================================================
[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);


% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
Hc  = spm_csd_fmri_mtf(Qp,DCM.M,DCM.U);                  % prediction
Ec  = spm_unvec(spm_vec(DCM.xY.y) - spm_vec(Hc),Hc);     % prediction error
 
 
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
M             = rmfield(DCM.M,'U'); 
M.dipfit.type = 'LFP';

M.U         = 1; 
M.l         = Ns;
qp          = Qp;
qp.L        = ones(1,Ns);             % set virtual electrode gain to unity
qp.b        = qp.b - 32;              % and suppress non-specific and
qp.c        = qp.c - 32;              % specific channel noise

[Hs Hz dtf] = spm_csd_mtf(qp,M,DCM.xU);
[ccf pst]   = spm_csd2ccf(Hs,DCM.M.Hz);
[coh fsd]   = spm_csd2coh(Hs,DCM.M.Hz);
DCM.dtf     = dtf;
DCM.ccf     = ccf;
DCM.coh     = coh;
DCM.fsd     = fsd;
DCM.pst     = pst;
DCM.Hz      = Hz;

 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM.Ep = Qp;                   % conditional expectation
DCM.Cp = Cp;                   % conditional covariance
DCM.Pp = Pp;                   % conditional probability
DCM.Hc = Hc;                   % conditional responses (y), channel space
DCM.Rc = Ec;                   % conditional residuals (y), channel space
DCM.Hs = Hs;                   % conditional responses (y), source space
DCM.Ce = exp(-Eh);             % ReML error covariance
DCM.F  = F;                    % Laplace log evidence
DCM.ID = ID;                   % data ID
 
 
save(DCM.name,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));

function DCM = spm_dcm_fmri_csd(P)
% Estimates parameters of a DCM using cross spectral fMRI densities
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
% DCM.Y.y                            % responses (over time)
% DCM.n                              % number of regions
% DCM.v                              % number of scans
%
% This routine estimates the (A and C) parameters of a dynamic causal model
% of fMRI responses, using the complex cross spectra under stationarity
% assumptions. The cross spectra are estimated from regional timeseries
% (the nodes of the DCM graph) using a Bayesian multivariate autoregressive
% model. The complex cross spectra are then fitted using linear systems
% theory in frequency space, under the simple assumption that the observed
% spectra are the predicted spectra plus some scale free fluctuations
% (noise). The characterisation of the model parameters can then be
% examined in terms of directed transfer functions, spectral density and
% crosscorrelation functions at the neuronal level - having accounted for
% variations in haemodynamics at each node.
%
% NB: if DCM.Y.y{i} is a cell array of multiple time series (e.g., sessions
% or subjects), this routine will use DCM.b  as constraints on the
% connectivity parameters that can change over sessions. The posterior
% estimates in DCM.Ep.B  then correspond to the session specific deviations
% from the average in DCM.Ep.A. The remaining results  pertain to the
% average connectivity.  This facility can be used to test for between
% session (or subject) effects with a subsequent application of parametric
% empirical Bayes (PEB), applied to the  field 'B'.
%
% see also: spm_dcm_estimate
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_fmri_csd.m 7497 2018-11-24 17:00:25Z karl $

SVNid = '$Rev: 7497 $';

% Load DCM structure
%--------------------------------------------------------------------------

if isstruct(P)
    DCM = P;
    try, DCM.name; catch, DCM.name = sprintf('DCM_%s',date); end
    P   = DCM.name;
else
    load(P)
end


% check options and specification
%==========================================================================
try, DCM.options.two_state;  catch, DCM.options.two_state  = 0;     end
try, DCM.options.stochastic; catch, DCM.options.stochastic = 0;     end
try, DCM.options.centre;     catch, DCM.options.centre     = 0;     end
try, DCM.options.analysis;   catch, DCM.options.analysis   = 'CSD'; end
try, DCM.options.order;      catch, DCM.options.order      = 8;     end
try, DCM.options.nograph;    catch, DCM.options.nograph    = spm('CmdLine');  end


% parameter initialisation
%--------------------------------------------------------------------------
try, DCM.M.P = DCM.options.P; end

% check max iterations
%--------------------------------------------------------------------------
try
    DCM.options.maxit;
catch    
    if isfield(DCM.options,'nN')
        DCM.options.maxit = DCM.options.nN;
        warning('options.nN is deprecated. Please use options.maxit');
    else
        DCM.options.maxit = 128;
    end
end

try DCM.M.Nmax; catch, DCM.M.Nmax = DCM.options.maxit; end

% check max nodes
%--------------------------------------------------------------------------
try
    DCM.options.maxnodes;
catch
    if isfield(DCM.options,'nmax')
        DCM.options.maxnodes = DCM.options.nmax;
        warning('options.nmax is deprecated. Please use options.maxnodes');
    else
        DCM.options.maxnodes = 8;
    end
end

% sizes
%--------------------------------------------------------------------------
try, DCM.U.u; catch, DCM.U.u = [];            end
try, DCM.n;   catch, DCM.n = size(DCM.a,1);   end
try, DCM.v;   catch, DCM.v = size(DCM.Y.y,1); end


% analysis and options
%--------------------------------------------------------------------------
DCM.options.induced    = 1;
DCM.options.nonlinear  = 0;
DCM.options.stochastic = 0;


% organise response variables: detrend outputs (and inputs)
%==========================================================================
DCM.Y.y = spm_detrend(DCM.Y.y);
if DCM.options.centre
    DCM.U.u = spm_detrend(DCM.U.u);
end

% scale timeseries to a precision of four
%--------------------------------------------------------------------------
vY          = spm_vec(DCM.Y.y);
scale       = 1/std(vY)/4;
DCM.Y.y     = spm_unvec(vY*scale,DCM.Y.y);
DCM.Y.scale = scale;

% disable high order parameters and check for models with no inputs
%--------------------------------------------------------------------------
n       = DCM.n;
if iscell(DCM.Y.y)
    
    % augment with between session constraints
    %----------------------------------------------------------------------
    if size(DCM.b,3) ~= numel(DCM.Y.y)
        for i = 1:numel(DCM.Y.y)
            DCM.b(:,:,i)  = DCM.b(:,:,1);
        end
    end
else
    DCM.b      = zeros(n,n,0);
end
if isempty(DCM.c) || isempty(DCM.U.u)
    DCM.c      = zeros(n,0);
    DCM.U.u    = zeros(DCM.v,1);
    DCM.U.name = {'null'};
end
DCM.d   = zeros(n,n,0);


% priors (and initial states)
%==========================================================================
[pE,pC,x] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);

% eigenvector constraints on pC for large models
%--------------------------------------------------------------------------
if n > DCM.options.maxnodes
    
    % remove confounds and find principal (maxnodes) modes
    %----------------------------------------------------------------------
    try
        y   = DCM.Y.y - DCM.Y.X0*(pinv(DCM.Y.X0)*DCM.Y.y);
    catch
        y   = spm_detrend(DCM.Y.y);
    end
    V       = spm_svd(y');
    V       = V(:,1:DCM.options.maxnodes);
    
    % remove minor modes from priors on A
    %----------------------------------------------------------------------
    j       = 1:(n*n);
    V       = kron(V*V',V*V');
    pC(j,j) = V*pC(j,j)*V';
    
end

% place eigenmodes in model if DCM.a is a vector (of eigenvalues)
%--------------------------------------------------------------------------
if isvector(DCM.a)
    DCM.M.modes = spm_svd(cov(DCM.Y.y));
end

% check for pre-specified priors
%--------------------------------------------------------------------------
hE       = 8;
hC       = 1/128;
try, pE  = DCM.M.pE; pC  = DCM.M.pC; end
try, hE  = DCM.M.hE; hC  = DCM.M.hC; end

% create DCM
%--------------------------------------------------------------------------
DCM.M.IS = 'spm_csd_fmri_mtf';
DCM.M.g  = @spm_gx_fmri;
DCM.M.f  = @spm_fx_fmri;
DCM.M.x  = x;
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.hE = hE;
DCM.M.hC = hC;
DCM.M.n  = length(spm_vec(x));
DCM.M.m  = size(DCM.U.u,2);
DCM.M.l  = n;
DCM.M.p  = DCM.options.order;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM.M.u  = sparse(n,1);

% get data-features (MAR(p) model)
%==========================================================================
DCM.Y.p  = DCM.M.p;
DCM      = spm_dcm_fmri_csd_data(DCM);
DCM.M.Hz = DCM.Y.Hz;
DCM.M.dt = 1/2;
DCM.M.N  = 32;
DCM.M.ns = 1/DCM.Y.dt;


% scale input (to a variance of 1/8)
%--------------------------------------------------------------------------
if any(diff(DCM.U.u))
    ccf         = spm_csd2ccf(DCM.U.csd,DCM.Y.Hz);
    DCM.U.scale = max(spm_vec(ccf))*8;
    DCM.U.csd   = spm_unvec(spm_vec(DCM.U.csd)/DCM.U.scale,(DCM.U.csd));
end


% complete model specification and invert
%==========================================================================

% precision of spectral observation noise
%--------------------------------------------------------------------------
DCM.Y.Q  = spm_dcm_csd_Q(DCM.Y.csd);
DCM.Y.X0 = sparse(size(DCM.Y.Q,1),0);
DCM.Y.p  = DCM.M.p;

% Variational Laplace: model inversion (using spectral responses)
%==========================================================================
Y.y          = DCM.Y.csd;
[Ep,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.U,Y);


% Bayesian inference {threshold = prior}
%--------------------------------------------------------------------------
warning('off','SPM:negativeVariance');
dp     = spm_vec(Ep) - spm_vec(pE);
Pp     = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Ep);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd) and error (haemodynamics)
%==========================================================================
Hc     = spm_csd_fmri_mtf(Ep,DCM.M,DCM.U);           % prediction
Ec     = spm_unvec(spm_vec(Y.y) - spm_vec(Hc),Hc);   % prediction error

% Bilinear representation and first-order hemodynamic kernel
%--------------------------------------------------------------------------
M      = DCM.M;                                      % model
Qp     = Ep;                                         % posterior parameters
Qp.B   = zeros(n,n,0);                               % no session effects
Qp.C   = speye(n,n);                                 % Switch to endogenous
[S,H1] = spm_dcm_mtf(Qp,M);                          % haemodynamic kernel

% and neuronal kernels
%==========================================================================
M.g    = @(x,u,P,M) x(:,1);                          % neuronal observer
[S,K1] = spm_dcm_mtf(Qp,M);                          % neuronal kernel

% predictions (at the level of neuronal interactions)
%--------------------------------------------------------------------------
Qp.b        = Qp.b - 32;                             % Switch off noise
Qp.c        = Qp.c - 32;                             % Switch off noise
Qp.C        = Ep.C;
[Hs,Hz,dtf] = spm_csd_fmri_mtf(Qp,M,DCM.U);
[ccf,pst]   = spm_csd2ccf(Hs,Hz);
[coh,fsd]   = spm_csd2coh(Hs,Hz);
DCM.dtf     = dtf;
DCM.ccf     = ccf;
DCM.coh     = coh;
DCM.fsd     = fsd;
DCM.pst     = pst;
DCM.Hz      = Hz;
DCM.H1      = H1;
DCM.K1      = K1;
 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM.Ep = Ep;                   % conditional expectation
DCM.Cp = Cp;                   % conditional covariance
DCM.Pp = Pp;                   % conditional probability
DCM.Hc = Hc;                   % conditional responses (y), BOLD space
DCM.Rc = Ec;                   % conditional residuals (y), BOLD space
DCM.Hs = Hs;                   % conditional responses (y), neural space
DCM.Ce = exp(-Eh);             % ReML error covariance
DCM.F  = F;                    % Laplace log evidence

% Save SPM version and revision number of code used
%--------------------------------------------------------------------------
[DCM.version.SPM.version, DCM.version.SPM.revision] = spm('Ver');
DCM.version.DCM.version  = spm_dcm_ui('Version');
DCM.version.DCM.revision = SVNid;

% and save
%--------------------------------------------------------------------------
save(P,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));

return












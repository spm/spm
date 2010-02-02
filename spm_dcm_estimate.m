function [DCM] = spm_dcm_estimate(P)
% Estimates parameters of a DCM (bilinear or nonlinear) for fMRI data
% FORMAT [DCM] = spm_dcm_estimate(DCM)
%   DCM - the DCM or its filename
%
% Expects
%--------------------------------------------------------------------------
% DCM.a;                             % switch on endogenous connections
% DCM.b;                             % switch on bilinear modulations
% DCM.c;                             % switch on exogenous connections
% DCM.U;                             % exogenous inputs
% DCM.Y;                             % responses
% DCM.Y.X0;                          % confounds
% DCM.n;                             % number of regions
% DCM.v;                             % number of scans
%
% Options
%--------------------------------------------------------------------------
% DCM.options.two_state              % two regional populations (E and I)
% DCM.options.stochastic             % fluctuations on hidden states
% DCM.options.nonlinear              % interactions among hidden states
%
% Evaluates:
%--------------------------------------------------------------------------
% DCM.M                              % Model structure
% DCM.Ep                             % Condition means (parameter structure)
% DCM.Cp                             % Conditional covariances
% DCM.Vp                             % Conditional variances
% DCM.Pp                             % Conditional probabilities
% DCM.H1                             % 1st order hemodynamic kernels
% DCM.H2                             % 2nd order hemodynamic kernels
% DCM.K1                             % 1st order neuronal kernels
% DCM.K2                             % 2nd order neuronal kernels
% DCM.R                              % residuals
% DCM.y                              % predicted data
% DCM.T                              % Threshold for Posterior inference
% DCM.Ce                             % Error variance for each region
% DCM.F                              % Free-energy bound on log evidence
% DCM.ID                             % Data ID
% DCM.AIC                            % Akaike Information criterion
% DCM.BIC                            % Bayesian Information criterion
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_estimate.m 3708 2010-02-02 20:15:31Z karl $
 
 
% load DCM structure
%--------------------------------------------------------------------------
if ~nargin
 
    %-display model details
    %----------------------------------------------------------------------
    Finter = spm_figure('GetWin','Interactive');
    set(Finter,'name','Dynamic Causal Modeling')
 
    %-get DCM
    %----------------------------------------------------------------------
    P = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    spm('Pointer','Watch')
    spm('FigName','Estimation in progress');
 
end
if isstruct(P)
    DCM = P;
    P   = ['DCM-' date];
else
    load(P)
end
 
% check options
%==========================================================================
try, DCM.options.two_state;  catch, DCM.options.two_state  = 0; end
try, DCM.options.stochastic; catch, DCM.options.stochastic = 0; end
try, DCM.options.nonlinear;  catch, DCM.options.nonlinear  = 0; end
 
% unpack
%--------------------------------------------------------------------------
U  = DCM.U;                             % exogenous inputs
Y  = DCM.Y;                             % responses
n  = DCM.n;                             % number of regions
v  = DCM.v;                             % number of scans
 
 
% check scaling of Y (enforcing a maximum change of 4%
%--------------------------------------------------------------------------
scale   = max(max((Y.y))) - min(min((Y.y)));
scale   = 4/max(scale,4);
Y.y     = Y.y*scale;
Y.scale = scale;

% normalise inputs
%--------------------------------------------------------------------------
U.u  = spm_detrend(U.u);
 
% check DCM.d (for nonlinear DCMs)
%--------------------------------------------------------------------------
try
    DCM.d;
catch
    DCM.d = zeros(n,n,0);
    DCM.options.nonlinear = 0;
end
 
% check confounds (add constant if necessary)
%--------------------------------------------------------------------------
if ~size(Y.X0,2), Y.X0 = ones(v,1); end
 
 
% specify parameters for spm_int_D (ensuring updates every second or so)
%--------------------------------------------------------------------------
if DCM.options.nonlinear
    M.IS     = 'spm_int_D';
    M.nsteps = round(max(Y.dt,1));
    M.states = 1:n;
else
    M.IS  = 'spm_int';
end
 
% priors
%--------------------------------------------------------------------------
if DCM.options.two_state
    [pE,pC] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,'2s');
    x       = sparse(n,6);
else
    [pE,pC] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d);
    x       = sparse(n,5);
end
 
% fMRI slice time sampling
%--------------------------------------------------------------------------
try, M.delays = DCM.delays; end
try, M.TE     = DCM.TE;     end
 
 
% complete model specification
%--------------------------------------------------------------------------
M.f   = 'spm_fx_fmri';
M.g   = 'spm_gx_fmri';
M.x   = x;
M.pE  = pE;
M.pC  = pC;
M.m   = size(U.u,2);
M.n   = size(x(:),1);
M.l   = size(x,1);
M.N   = 32;
M.dt  = 16/M.N;
M.ns  = v;
 
 
% nonlinear system identification (nlsi)
%==========================================================================
 
% nonlinear system identification (Variational EM) - deterministic DCM
%--------------------------------------------------------------------------
[Ep,Cp,Eh,F]  = spm_nlsi_GN(M,U,Y);
 
% proceed to stochastic (initialising with determinidstic estimates)
%--------------------------------------------------------------------------
if DCM.options.stochastic
    
    
    % specify bilinear approximation scheme for spm_DEM_eval
    %----------------------------------------------------------------------
    if size(pE.A,3), M.E.linear = 1; end   % linear model
    if size(pE.B,3), M.E.linear = 2; end   % bilinear model
    if size(pE.D,3), M.E.linear = 3; end   % nonlinear models
 
    % Decimate U.u from micro-time
    % ---------------------------------------------------------------------
    u       = U.u;
    y       = Y.y;
    Dy      = spm_dctmtx(size(y,1),size(y,1));
    Du      = spm_dctmtx(size(u,1),size(y,1));
    Dy      = Dy*sqrt(size(y,1)/size(u,1));
    u       = Dy*(Du'*u);
 
    % DEM Structure: place model, data, input and confounds in DEM
    % ---------------------------------------------------------------------
    DEM.M   = M;
    DEM.Y   = y';
    DEM.U   = u';
    DEM.X   = Y.X0';
 
    % set inversion parameters
    % ---------------------------------------------------------------------
    DEM.M(1).E.s  = 1/8;        % smoothness of random fluctuations
    DEM.M(1).E.d  = 2;          % embedding dimension 
    DEM.M(1).E.n  = 4;          % embedding dimension
    DEM.M(1).E.nN = 8;          % maximum number of DEM iterations
 
    % adjust M.f (DEM works in time bins not seconds) and initialize M.P
    % ---------------------------------------------------------------------
    DEM.M(1).f  = inline([M.f '(x,v,P)*' num2str(Y.dt)],'x','v','P');
    DEM.M(1).P  = Ep;
    
    % delays (DEM works in time bins not seconds)
    % ---------------------------------------------------------------------
    DEM.M(1).delays = M.delays/Y.dt;
    
    % Precision on hidden states
    % ---------------------------------------------------------------------
    DEM.M(1).xP = 32;
 
    
    % Specify hyper-priors on precisions: level 1
    % ---------------------------------------------------------------------
    DEM.M(1).Q  = spm_Ce(ones(1,n));
    DEM.M(1).hE = Eh;          % prior expectation of log precision (noise)
    DEM.M(1).hC = 8;           % prior covariance  of log precision (noise)
    DEM.M(1).gE = 8;           % prior expectation of log precision (state)
    DEM.M(1).gC = 1;           % prior covariance  of log precision (state)
 
    % and level 2
    %----------------------------------------------------------------------
    DEM.M(2).V  = 4;           % prior expectation of log precision (cause)
 
 
    % Generalised filtering (under the Laplace assumption)
    % =====================================================================
    DEM    = spm_LAP(DEM);
 
    % Save DEM estimates
    %----------------------------------------------------------------------
    DCM.qU = DEM.qU;
    DCM.qP = DEM.qP;
    DCM.qH = DEM.qH;
    DCM.Fd = F;
    
    % unpack results
    % ---------------------------------------------------------------------
    F      = DEM.F(end);
    Ep     = DEM.qP.P{1};
    Cp     = DEM.qP.C;
 
    % predicted responses (y) and residuals (R) from DEM
    %--------------------------------------------------------------------------
    y      = DEM.qU.v{1}';
    R      = DEM.qU.z{1}';
    R      = R - Y.X0*inv(Y.X0'*Y.X0)*(Y.X0'*R);
    Ce     = exp(-DEM.qH.h{1});
 
else
 
    % predicted responses (y) and residuals (R) from EM
    %--------------------------------------------------------------------------
    y      = feval(M.IS,Ep,M,U);
    R      = Y.y - y;
    R      = R - Y.X0*inv(Y.X0'*Y.X0)*(Y.X0'*R);
    Ce     = exp(-Eh);
    
end
 
 
% Bilinear representation and first-order hemodynamic kernel
%--------------------------------------------------------------------------
[M0,M1,L1,L2] = spm_bireduce(M,Ep);
[H0,H1] = spm_kernels(M0,M1,L1,L2,M.N,M.dt);
 
% and neuronal kernels
%--------------------------------------------------------------------------
L       = sparse(1:n,[1:n] + 1,1,n,length(M0));
[K0,K1] = spm_kernels(M0,M1,L,M.N,M.dt);
 
 
% Bayesian inference and variance {threshold T = 0}
%--------------------------------------------------------------------------
T       = 0;
sw      = warning('off','SPM:negativeVariance');
Pp      = spm_unvec(1 - spm_Ncdf(T,abs(spm_vec(Ep)),diag(Cp)),Ep);
Vp      = spm_unvec(diag(Cp),Ep);
warning(sw);
 
 
% Store parameter estimates
%--------------------------------------------------------------------------
DCM.M   = M;
DCM.Y   = Y;
DCM.U   = U;
DCM.Ce  = Ce;
DCM.Ep  = Ep;
DCM.Cp  = Cp;
DCM.Pp  = Pp;
DCM.Vp  = Vp;
DCM.H1  = H1;
DCM.K1  = K1;
DCM.R   = R;
DCM.y   = y;
DCM.T   = T;
 
 
% Data ID and log-evidence
%==========================================================================
if isfield(M,'FS')
    try
        ID = spm_data_id(feval(M.FS,Y.y,M));
    catch
        ID = spm_data_id(feval(M.FS,Y.y));
    end
else
    ID     = spm_data_id(Y.y);
end
 
% Save approximations to model evidence: negative free energy, AIC, BIC
%--------------------------------------------------------------------------
evidence   = spm_dcm_evidence(DCM);
DCM.F      = F;
DCM.ID     = ID;
DCM.AIC    = evidence.aic_overall;
DCM.BIC    = evidence.bic_overall;
 
%-Save DCM
%--------------------------------------------------------------------------
save(P,'DCM','F','Ep','Cp');
 
if ~nargin
    spm('Pointer','Arrow');
    spm('FigName','Done');
end

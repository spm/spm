function [DCM] = spm_dcm_estimate(P)
% Estimates parameters of a DCM (bilinear or nonlinear) for fMRI data
% FORMAT [DCM] = spm_dcm_estimate(DCM)
%   DCM - DCM structure or its filename
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
%
% Options
%--------------------------------------------------------------------------
% DCM.options.two_state              % two regional populations (E and I)
% DCM.options.stochastic             % fluctuations on hidden states
% DCM.options.centre                 % mean-centre inputs
% DCM.options.nonlinear              % interactions among hidden states
% DCM.options.nograph                % graphical display
% DCM.options.induced                % switch for CSD data features
% DCM.options.P                      % starting estimates for parameters
% DCM.options.hidden                 % indices of hidden regions
% DCM.options.maxnodes               % maximum number of (effective) nodes
% DCM.options.maxit                  % maximum number of iterations
% DCM.options.hE                     % expected precision of the noise
% DCM.options.hC                     % variance of noise expectation
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
% Copyright (C) 2002-2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_estimate.m 7479 2018-11-09 14:17:33Z peter $
 
 
SVNid = '$Rev: 7479 $';
 
%-Load DCM structure
%--------------------------------------------------------------------------
if ~nargin
    
    %-Display model details
    %----------------------------------------------------------------------
    Finter = spm_figure('GetWin','Interactive');
    set(Finter,'name','Dynamic Causal Modelling')
    
    %-Get DCM
    %----------------------------------------------------------------------
    [P, sts] = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, DCM = []; return; end
    spm('Pointer','Watch')
    spm('FigName','Estimation in progress');
    
end
 
if isstruct(P)
    DCM = P;
else
    load(P)
end
 
% check options
%==========================================================================
try, DCM.options.two_state;  catch, DCM.options.two_state  = 0;     end
try, DCM.options.stochastic; catch, DCM.options.stochastic = 0;     end
try, DCM.options.nonlinear;  catch, DCM.options.nonlinear  = 0;     end
try, DCM.options.centre;     catch, DCM.options.centre     = 0;     end
try, DCM.options.hidden;     catch, DCM.options.hidden     = [];    end
try, DCM.options.hE;         catch, DCM.options.hE         = 6;     end
try, DCM.options.hC;         catch, DCM.options.hC         = 1/128; end
try, DCM.n;                  catch, DCM.n = size(DCM.a,1);          end
try, DCM.v;                  catch, DCM.v = size(DCM.Y.y,1);        end
 
try, M.nograph = DCM.options.nograph; catch, M.nograph = spm('CmdLine');end
 
% check max iterations
%--------------------------------------------------------------------------
try
    DCM.options.maxit;
catch    
    if isfield(DCM.options,'nN')
        DCM.options.maxit = DCM.options.nN;
        warning('options.nN is deprecated. Please use options.maxit');
    elseif DCM.options.stochastic
        DCM.options.maxit = 32;
    else
        DCM.options.maxit = 128;
    end
end
 
try M.Nmax = DCM.M.Nmax; catch, M.Nmax = DCM.options.maxit; end
 
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
 
% analysis and options
%--------------------------------------------------------------------------
DCM.options.induced  = 0;
 
% unpack
%--------------------------------------------------------------------------
U  = DCM.U;                             % exogenous inputs
Y  = DCM.Y;                             % responses
n  = DCM.n;                             % number of regions
v  = DCM.v;                             % number of scans
 
% detrend outputs (and inputs)
%--------------------------------------------------------------------------
Y.y = spm_detrend(Y.y);
if DCM.options.centre
    U.u = spm_detrend(U.u);
end
 
% check scaling of Y (enforcing a maximum change of 4%
%--------------------------------------------------------------------------
scale   = max(max((Y.y))) - min(min((Y.y)));
scale   = 4/max(scale,4);
Y.y     = Y.y*scale;
Y.scale = scale;
 
% check confounds (add constant if necessary)
%--------------------------------------------------------------------------
if ~isfield(Y,'X0'),Y.X0 = ones(v,1); end
if ~size(Y.X0,2),   Y.X0 = ones(v,1); end
 
% fMRI slice time sampling
%--------------------------------------------------------------------------
try, M.delays = DCM.delays; catch, M.delays = ones(n,1); end
try, M.TE     = DCM.TE;     end
 
% create priors
%==========================================================================
 
% check DCM.d (for nonlinear DCMs)
%--------------------------------------------------------------------------
try
    DCM.options.nonlinear = logical(size(DCM.d,3));
catch
    DCM.d = zeros(n,n,0);
    DCM.options.nonlinear = 0;
end
 
% specify parameters for spm_int_D (ensuring updates every second or so)
%--------------------------------------------------------------------------
if DCM.options.nonlinear
    M.IS     = 'spm_int_D';
    M.nsteps = round(max(Y.dt,1));
    M.states = 1:n;
else
    M.IS     = 'spm_int';
end
 
% check for endogenous DCMs, with no exogenous driving effects
%--------------------------------------------------------------------------
if isempty(DCM.c) || isempty(U.u)
    DCM.c  = zeros(n,1);
    DCM.b  = zeros(n,n,1);
    U.u    = zeros(v,1);
    U.name = {'null'};
end
if ~any(spm_vec(U.u)) || ~any(spm_vec(DCM.c))
    DCM.options.stochastic = 1;
end
 
% priors (and initial states)
%--------------------------------------------------------------------------
[pE,pC,x]  = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);
str        = 'Using specified priors ';
str        = [str '(any changes to DCM.a,b,c,d will be ignored)\n'];
 
try, M.P   = DCM.options.P;                end      % initial parameters
try, pE    = DCM.options.pE; fprintf(str); end      % prior expectation
try, pC    = DCM.options.pC; fprintf(str); end      % prior covariance
 
try, M.P   = DCM.M.P;                end            % initial parameters
try, pE    = DCM.M.pE; fprintf(str); end            % prior expectation
try, pC    = DCM.M.pC; fprintf(str); end            % prior covariance
 
 
% eigenvector constraints on pC for large models
%--------------------------------------------------------------------------
if n > DCM.options.maxnodes
    
    % remove confounds and find principal (nmax) modes
    %----------------------------------------------------------------------
    y       = Y.y - Y.X0*(pinv(Y.X0)*Y.y);
    V       = spm_svd(y');
    V       = V(:,1:DCM.options.maxnodes);
    
    % remove minor modes from priors on A
    %----------------------------------------------------------------------
    j       = 1:(n*n);
    V       = kron(V*V',V*V');
    pC(j,j) = V*pC(j,j)*V';
    
end
 
% hyperpriors over precision - expectation and covariance
%--------------------------------------------------------------------------
hE      = sparse(n,1) + DCM.options.hE;
hC      = speye(n,n)  * DCM.options.hC;
i       = DCM.options.hidden;
hE(i)   = -4;
hC(i,i) = exp(-16);
 
% complete model specification
%--------------------------------------------------------------------------
M.f  = 'spm_fx_fmri';                     % equations of motion
M.g  = 'spm_gx_fmri';                     % observation equation
M.x  = x;                                 % initial condition (states)
M.pE = pE;                                % prior expectation (parameters)
M.pC = pC;                                % prior covariance  (parameters)
M.hE = hE;                                % prior expectation (precisions)
M.hC = hC;                                % prior covariance  (precisions)
M.m  = size(U.u,2);
M.n  = size(x(:),1);
M.l  = size(x,1);
M.N  = 64;
M.dt = 32/M.N;
M.ns = v;

% nonlinear system identification (nlsi)
%==========================================================================
if ~DCM.options.stochastic
    
    % nonlinear system identification (Variational EM) - deterministic DCM
    %----------------------------------------------------------------------
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);
    
    % predicted responses (y) and residuals (R)
    %----------------------------------------------------------------------
    y      = feval(M.IS,Ep,M,U);
    R      = Y.y - y;
    R      = R - Y.X0*spm_inv(Y.X0'*Y.X0)*(Y.X0'*R);
    Ce     = exp(-Eh);   
    
else
    
    % proceed to stochastic (initialising with deterministic estimates)
    %======================================================================
    
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
    DEM.M(1).E.form = 'Gaussian';         % form of random fluctuations
    DEM.M(1).E.s    = 1/2;                % smoothness of fluctuations
    DEM.M(1).E.d    = 2;                  % embedding dimension
    DEM.M(1).E.n    = 4;                  % embedding dimension
    DEM.M(1).E.nN   = DCM.options.maxit;  % maximum number of iterations
    
    % adjust M.f (DEM works in time bins not seconds) and initialize M.P
    % ---------------------------------------------------------------------
    DEM.M(1).delays = M.delays/Y.dt;
    DEM.M(1).f      = inline([M.f '(x,v,P)*' num2str(Y.dt)],'x','v','P');
  
    
    % Specify hyper-priors on (log-precision of) observation noise
    % ---------------------------------------------------------------------
    DEM.M(1).Q  = spm_Ce(ones(1,n));      % precision components
    DEM.M(1).hE = hE;                     % prior expectation
    DEM.M(1).hC = hC;                     % prior covariance
    
    
    % allow (only) neuronal [x, s, f, q, v] hidden states to fluctuate
    % ---------------------------------------------------------------------
    W           = ones(n,1)*[12 16 16 16 16];
    DEM.M(1).xP = exp(6);                 % precision (hidden-state)
    DEM.M(1).W  = diag(exp(W));           % precision (hidden-motion)
    DEM.M(2).V  = exp(16);                % precision (hidden-cause)
    
    
    % Generalised filtering (under the Laplace assumption)
    % =====================================================================
    DEM    = spm_LAP(DEM);
    
    
    % Save DEM estimates
    %----------------------------------------------------------------------
    DCM.qU = DEM.qU;
    DCM.qP = DEM.qP;
    DCM.qH = DEM.qH;
    
    % unpack results
    % ---------------------------------------------------------------------
    F      = DEM.F(end);
    Ep     = DEM.qP.P{1};
    Cp     = DEM.qP.C;
    
    % predicted responses (y) and residuals (R)
    %----------------------------------------------------------------------
    y      = DEM.qU.v{1}';
    R      = DEM.qU.z{1}';
    R      = R - Y.X0*spm_inv(Y.X0'*Y.X0)*(Y.X0'*R);
    Ce     = exp(-DEM.qH.h{1});
        
end
 
 
% Bilinear representation and first-order hemodynamic kernel
%--------------------------------------------------------------------------
[M0,M1,L1,L2] = spm_bireduce(M,Ep);
[H0,H1] = spm_kernels(M0,M1,L1,L2,M.N,M.dt);
 
% and neuronal kernels
%--------------------------------------------------------------------------
L       = sparse(1:n,(1:n) + 1,1,n,length(M0));
[K0,K1] = spm_kernels(M0,M1,L,M.N,M.dt);
 
 
% Bayesian inference and variance {threshold: prior mean plus T = 0}
%--------------------------------------------------------------------------
T       = full(spm_vec(pE));
sw      = warning('off','SPM:negativeVariance');
Pp      = spm_unvec(1 - spm_Ncdf(T,abs(spm_vec(Ep)),diag(Cp)),Ep);
Vp      = spm_unvec(full(diag(Cp)),Ep);
warning(sw);
 
try,  M = rmfield(M,'nograph'); end
 
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
DCM.T   = 0;
 
 
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
 
% Save SPM version and revision number of code used
%--------------------------------------------------------------------------
[DCM.version.SPM.version, DCM.version.SPM.revision] = spm('Ver');
DCM.version.DCM.version  = spm_dcm_ui('Version');
DCM.version.DCM.revision = SVNid;
 
%-Save DCM
%--------------------------------------------------------------------------
if ~isstruct(P)
    save(P,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));
end
 
if ~nargin
    spm('Pointer','Arrow');
    spm('FigName','Done');
end

function [DCM] = spm_dcm_fmri_nmm(P)
% Estimates parameters of a DCM (neural mass model) for fMRI data
% FORMAT [DCM] = spm_dcm_estimate(DCM)
%   DCM - DCM structure or its filename
%
% Expects
%--------------------------------------------------------------------------
% DCM.a                              % switch on endogenous connections
% DCM.b                              % switch on bilinear modulations
% DCM.c                              % switch on exogenous connections
% DCM.U                              % exogenous inputs
% DCM.Y.y                            % responses
% DCM.Y.X0                           % confounds
% DCM.Y.Q                            % array of precision components
% DCM.n                              % number of regions
% DCM.v                              % number of scans
%
% Options
%--------------------------------------------------------------------------
% DCM.options.nmm                    % neural mass model
% DCM.options.centre                 % mean-centre inputs
% DCM.options.nograph                % graphical display
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
% $Id: spm_dcm_fmri_nmm.m 6856 2016-08-10 17:55:05Z karl $
 
 
SVNid = '$Rev: 6856 $';
 
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
    P   = ['DCM-' date '.mat'];
else
    load(P)
end
 
% check options
%==========================================================================
try, DCM.options.nmm;        catch, DCM.options.nmm        = 'TFM'; end
try, DCM.options.two_state;  catch, DCM.options.two_state  = 0;     end
try, DCM.options.stochastic; catch, DCM.options.stochastic = 0;     end
try, DCM.options.nonlinear;  catch, DCM.options.nonlinear  = 0;     end
try, DCM.options.centre;     catch, DCM.options.centre     = 0;     end
try, DCM.options.hidden;     catch, DCM.options.hidden     = [];    end
try, DCM.options.hE;         catch, DCM.options.hE         = 6;     end
try, DCM.options.hC;         catch, DCM.options.hC         = 1/128; end
try, DCM.options.maxit;      catch, DCM.options.maxit      = 32;    end
try, DCM.options.maxnodes;   catch, DCM.options.maxnodes   = 16;    end

try, DCM.n;                  catch, DCM.n = size(DCM.c,1);          end
try, DCM.v;                  catch, DCM.v = size(DCM.Y.y,1);        end
 
try, M.nograph = DCM.options.nograph; catch, M.nograph = spm('CmdLine'); end 
try, M.Nmax    = DCM.M.Nmax;          catch, M.Nmax = DCM.options.maxit; end
 
 
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
 
 
% priors (and initial states) - neuronal
%--------------------------------------------------------------------------
[nE,nC] = spm_dcm_neural_priors(DCM.a,DCM.b,DCM.c,DCM.options.nmm);
[xn,fn] = spm_dcm_x_neural(nE,DCM.options.nmm);

% fix some neuronal (neural mass) parameters
%--------------------------------------------------------------------------
p = {'T','D','R','G'};
for i = 1:numel(p)
    try
        nC.(p{i}) = spm_zeros(nC.(p{i}));
    end
end

pE.N = rmfield(nE,{'M','N','E','F'});
pC.N = rmfield(nC,{'M','N','E','F'});

% add neurovascular parameters
%--------------------------------------------------------------------------
pE.J       = sparse(4,3);   pC.J      = sparse(4,3) + 1/16;

% and hemodynamic priors
%--------------------------------------------------------------------------
bE.transit = sparse(n,1);  bC.transit = sparse(n,1) + 1/256;
bE.decay   = sparse(n,1);  bC.decay   = sparse(n,1) + 1/256;
bE.epsilon = sparse(1,1);  bC.epsilon = sparse(1,1) + 1/256;
x.x        = xn;
x.h        = sparse(n,4);

pE.H = bE;
pC.H = bC;


% Using specified priors
%--------------------------------------------------------------------------
str      = 'Using specified priors';
try, M.P = DCM.options.P;                end      % initial parameters
try, pE  = DCM.options.pE; fprintf(str); end      % prior expectation
try, pC  = DCM.options.pC; fprintf(str); end      % prior covariance



% hyperpriors over precision - expectation and covariance
%--------------------------------------------------------------------------
hE       = sparse(n,1) + DCM.options.hE;
hC       = speye(n,n)  * DCM.options.hC;
i        = DCM.options.hidden;
hE(i)    = -4;
hC(i,i)  = exp(-16);
 
% complete model specification
%--------------------------------------------------------------------------
M.IS = @spm_gen_fmri;                     % integration scheme
M.fn = fn;                                % neuronal equations of motion
M.x  = x;                                 % initial condition (states)
M.pE = pE;                                % prior expectation (parameters)
M.pC = pC;                                % prior covariance  (parameters)
M.hE = hE;                                % prior expectation (precisions)
M.hC = hC;                                % prior covariance  (precisions)
M.m  = size(U.u,2);
M.n  = size(spm_vec(x),1);
M.l  = n;
M.ns = v;

    
% nonlinear system identification (Variaitonal Laplace)
%==========================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);

% predicted responses (y) and residuals (R)
%--------------------------------------------------------------------------
y     = feval(M.IS,Ep,M,U);
R     = Y.y - y;
R     = R - Y.X0*spm_inv(Y.X0'*Y.X0)*(Y.X0'*R);
Ce    = exp(-Eh);
    
% Bilinear representation and first-order hemodynamic kernel
%--------------------------------------------------------------------------
H.f   = @spm_fx_hdm;
H.g   = @spm_gx_hdm;
H.x   = M.x.h;
H.m   = M.m;

[H0,H1] = spm_kernels(H,Ep.H,64,1/2);
 
% and neuronal kernels
%--------------------------------------------------------------------------
N.f   = M.fn;
N.x   = M.x.x;
N.m   = M.m;

[K0,K1] = spm_kernels(N,Ep.N,64,8/1000);
 
 
% Bayesian inference and variance {threshold: prior mean plus T = 0}
%--------------------------------------------------------------------------
T     = full(spm_vec(pE));
sw    = warning('off','SPM:negativeVariance');
Pp    = spm_unvec(1 - spm_Ncdf(T,abs(spm_vec(Ep)),diag(Cp)),Ep);
Vp    = spm_unvec(full(diag(Cp)),Ep);
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


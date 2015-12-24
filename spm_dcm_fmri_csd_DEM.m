function DCM = spm_dcm_fmri_csd_DEM(DCM)
% Estimates parameters of a DCM using cross spectral fMRI densities
% FORMAT DCM = spm_dcm_fmri_csd_DEM(DCM)
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
% This routine estimates the parameters of a hierarchical model
% of fMRI responses, using the complex cross spectra under stationarity
% assumptions. The cross spectra are estimated from regional timeseries
% (the nodes of the DCM graph) using a Bayesian multivariate autoregressive
% model. The complex cross spectra are then fitted using linear systems
% theory in frequency space, under the simple assumption that the observed
% spectra are the predicted spectra plus some smooth Gaussian fluctuations
% (noise). The characterisation of the model parameters can then be
% examined in terms of directed transfer functions, spectral density and
% crosscorrelation functions at the neuronal level - having accounted for
% variations in haemodynamics at each node.
%
% This scheming uses a hierarchical generative model of connectivity with
% hierarchical constraints on the edges and therefore uses the expectation
% and maximisation stepits of dynamic expectation maximisation. Here, the
% hidden causes at the first level are the effective connectivity and the
% hidden causes at the second level are the Lyapunov exponents or 
% eigenvalues of a symmetrical Jacobian or effective connectivity matrix:
% see DEM_demo_modes_fMRI.m
%
% see also: spm_dcm_estimate
%           spm_dcm_fmri_csd
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_csd_DEM.m 6656 2015-12-24 16:49:52Z guillaume $


% get DCM
%--------------------------------------------------------------------------
if ischar(DCM), load(DCM); end

% DCM name
%--------------------------------------------------------------------------
try, DCM.name; catch, DCM.name = sprintf('DCM_%s',date); end

% check options and specification
%==========================================================================
try, DCM.options.two_state;  catch, DCM.options.two_state  = 0;     end
try, DCM.options.stochastic; catch, DCM.options.stochastic = 0;     end
try, DCM.options.centre;     catch, DCM.options.centre     = 0;     end
try, DCM.options.nmax;       catch, DCM.options.nmax       = 8;     end
try, DCM.options.embedding;  catch, DCM.options.embedding  = 3;     end
try, DCM.options.nograph;    catch, DCM.options.nograph    = spm('CmdLine');  end


% parameter initialisation
%--------------------------------------------------------------------------
try, DCM.M.P = DCM.options.P;  end


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

% check scaling of Y (enforcing a maximum change of 4%
%--------------------------------------------------------------------------
scale       = max(max((DCM.Y.y))) - min(min((DCM.Y.y)));
scale       = 4/max(scale,4);
DCM.Y.y     = DCM.Y.y*scale;
DCM.Y.scale = scale;

% disable high order parameters and check for models with no inputs
%--------------------------------------------------------------------------
n       = DCM.n;
DCM.b   = zeros(n,n,0);
DCM.d   = zeros(n,n,0);
if isempty(DCM.c) || isempty(DCM.U.u)
    DCM.c      = zeros(DCM.n,1);
    DCM.U.u    = zeros(DCM.v,1);
    DCM.U.name = {'null'};
end

% priors (and initial states)
%==========================================================================
[pE,pC,x] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);


% eigenvector constraints on pC for large models
%--------------------------------------------------------------------------
if n > DCM.options.nmax
    
    % remove confounds and find principal (nmax) modes
    %----------------------------------------------------------------------
    try
        y   = DCM.Y.y - DCM.Y.X0*(pinv(DCM.Y.X0)*DCM.Y.y);
    catch
        y   = spm_detrend(DCM.Y.y);
    end
    V       = spm_svd(y');
    V       = V(:,1:DCM.options.nmax);
    
    % remove minor modes from priors on A
    %----------------------------------------------------------------------
    j       = 1:(n*n);
    V       = kron(V*V',V*V');
    pC(j,j) = V*pC(j,j)*V';
    
end


% create DCM
%--------------------------------------------------------------------------
DCM.M.IS = 'spm_csd_fmri_mtf';
DCM.M.FS = 'spm_fs_fmri_csd';
DCM.M.g  = @spm_gx_fmri;
DCM.M.f  = @spm_fx_fmri;
DCM.M.x  = x;
% DCM.M.pE = pE;
% DCM.M.pC = pC;
% DCM.M.hE = 6;
% DCM.M.hC = 1/128;
DCM.M.n  = length(spm_vec(x));
DCM.M.m  = size(DCM.U.u,2);
DCM.M.l  = n;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM.M.u  = sparse(n,1);

% get data-features (MAR(8) model)
%==========================================================================
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

% precision of spectral observation noise: AR(1/2)
%--------------------------------------------------------------------------
y        = feval(DCM.M.FS,DCM.Y.csd,DCM.M);
m        = size(y,2)*size(y,3);
q        = spm_Q(1/2,size(y,1),1);
Q        = kron(speye(m,m),q);
DCM.Y.Q  = Q;
DCM.Y.X0 = sparse(size(Q,1),0);


% global DCM
%==========================================================================
global GLOBAL_DCM
GLOBAL_DCM   = DCM;


% data (and priors)
%--------------------------------------------------------------------------
DEM.Y        = spm_vec(spm_fs_fmri_csd(DCM.Y.csd,DCM.M));
DEM.M.E.nD   = 8;
DEM.M.E.nE   = 8;

% functional modes and eigenvalues
%--------------------------------------------------------------------------
if 1
    
    % hidden causes as Lyapunov exponents or eigenvalues
    %----------------------------------------------------------------------
    P        = pE;
    P.modes  = spm_svd(cov(DCM.Y.y));
    try
        v    = DCM.options.v(1:DCM.options.embedding);
    catch
        v    = zeros(DCM.options.embedding,1);
    end
    
else
    
    % hidden causes as scaled eigenvectors
    %----------------------------------------------------------------------
    P        = pE;
    try
        v    = DCM.options.v(1:DCM.options.embedding,:);
    catch
        v    = spm_svd(corr(DCM.Y.y));
        v    = v(:,1:DCM.options.embedding)'/64;
    end
    
end

% generative model - for DEM
%==========================================================================
DEM.M(1).g   = 'spm_dcm_fmri_csd_gen';
DEM.M(1).hE  = 6;
DEM.M(1).hC  = exp(-6);
DEM.M(1).Q   = Q;
DEM.M(1).m   = length(pC);
DEM.M(1).n   = 0;
DEM.M(1).l   = length(Q);

DEM.M(2).g   = 'spm_dcm_fmri_graph_gen';
DEM.M(2).V   = spm_inv(pC);
DEM.M(2).v   = pE;
DEM.M(2).pE  = P;

DEM.M(3).V   = exp(0);
DEM.M(3).v   = v;

% Variational Laplace
%--------------------------------------------------------------------------
DEM    = spm_DEM(DEM);

% convert back into DCM format: [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);
%--------------------------------------------------------------------------
Ep     = DEM.qU.v{2};
Cp     = DEM.qU.C{1}(1:length(Ep),1:length(Ep));
Eh     = DEM.qH.h{1};
F      = DEM.F(end);
Ep     = spm_unvec(Ep,pE);

% Bayesian inference {threshold = prior}
%--------------------------------------------------------------------------
warning('off','SPM:negativeVariance');
dp     = spm_vec(Ep) - spm_vec(pE);
Pp     = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Ep);
warning('on', 'SPM:negativeVariance');


% predictions (csd) and error (haemodynamics)
%==========================================================================
Hc     = feval(DCM.M.IS,Ep,DCM.M,DCM.U);             % prediction
Ec     = DCM.Y.csd - Hc;                             % error

% Bilinear representation and first-order hemodynamic kernel
%--------------------------------------------------------------------------
M      = DCM.M;                                      % model
Qp     = Ep;                                         % posterior parameters
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
DCM.DEM     = DEM;

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

% and save
%--------------------------------------------------------------------------
save(DCM.name,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));

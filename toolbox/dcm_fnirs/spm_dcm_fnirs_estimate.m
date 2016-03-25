function [DCM] = spm_dcm_fnirs_estimate(P)
% Estimate parameters of a DCM for fNIRS data
% FORMAT [DCM] = spm_dcm_fnirs_estimate(P)
%
% P - Name of DCM file 
%
% DCM - DCM structure
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
% DCM.options.nmax                   % maximum number of (effective) nodes
% DCM.options.nN                     % maximum number of iterations
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
% Note: This code 
% (i) shows best results with spm_nlsi_GN.m (version 6481):
% Karl Friston $Id: spm_dcm_fnirs_estimate.m 6754 2016-03-25 06:44:58Z will $
% (ii) is based on spm_dcm_estimate.m:
% Karl Friston $Id: spm_dcm_fnirs_estimate.m 6754 2016-03-25 06:44:58Z will $
%__________________________________________________________________________
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_dcm_fnirs_estimate.m 6754 2016-03-25 06:44:58Z will $

%--------------------------------------------------------------------------
%-Load DCM structure
if ~nargin
 
    %----------------------------------------------------------------------
    %-Display model details
    Finter = spm_figure('GetWin','Interactive');
    set(Finter,'name','Dynamic Causal Modelling')
 
    %----------------------------------------------------------------------
    %-Get DCM
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

%--------------------------------------------------------------------------
% check options
try, DCM.options.two_state;  catch, DCM.options.two_state  = 0;     end
try, DCM.options.stochastic; catch, DCM.options.stochastic = 0;     end
try, DCM.options.nonlinear;  catch, DCM.options.nonlinear  = 0;     end
try, DCM.options.centre;     catch, DCM.options.centre     = 0;     end
try, DCM.options.nmax;       catch, DCM.options.nmax       = 8;     end
try, DCM.options.nN;         catch, DCM.options.nN         = 32;    end
try, DCM.options.hidden;     catch, DCM.options.hidden     = [];    end
try, DCM.n;                  catch, DCM.n = size(DCM.a,1);          end
try, DCM.v;                  catch, DCM.v = size(DCM.Y.y,1);        end

try, M.nograph = DCM.options.nograph; catch, M.nograph = spm('CmdLine');end

%--------------------------------------------------------------------------
% Unpack
U  = DCM.U;        % inputs
Y  = DCM.Y;       % responses
v  = DCM.v;        % number of scans
n = DCM.n;        % number of sources of interest 
m = size(U.u,2);  % number of inputs

%--------------------------------------------------------------------------
% load fNIRS measurements
ny = size(Y.P.fname.nirs, 1); 
nch = size(Y.P.rois, 2); 
for i = 1:ny
    raw = load(deblank(Y.P.fname.nirs{i, 1})); 
    y{i,1} = reshape(raw.Y.od(:, Y.P.rois, :), raw.P.ns, nch * size(Y.P.wav, 2)); 
end
Y.P.nch = nch; 

%--------------------------------------------------------------------------
% apply temporal filtering 
if strcmpi(Y.K.type, 'Butterworth IIR')
    addpath(fullfile(spm('Dir'), 'external', 'fieldtrip', 'preproc'));
    forder = 5; ftype = 'but'; fdir = 'twopass';
    
    nf = size(Y.K.cutoff, 1); 
    for i = 1:ny
        for j = 1:nf
            if Y.K.cutoff(j, 1) == 0 % high pass filter
                y{i,1} = ft_preproc_highpassfilter(y{i,1}', Y.P.fs, Y.K.cutoff(j, 2), forder, ftype, fdir)';
            else % band-stop filter
                y{i,1} = ft_preproc_bandstopfilter(y{i,1}', Y.P.fs , Y.K.cutoff(j, :), forder, ftype, fdir)';
            end
        end
    end
end

%--------------------------------------------------------------------------
% concatenate the time series 
y = cell2mat(y); 

%--------------------------------------------------------------------------
% average time series over trials 
if strcmpi(Y.W.type, 'y')
    nc = size(Y.W.onsets, 2); 
    for i = 1:nc, wy{i,1} = spm_fnirs_wavg(y, round(Y.W.onsets{i}.*Y.P.fs)+1, round(Y.W.durations(i).*Y.P.fs)); end 
    y = cell2mat(wy); 
end 
Y.P.ns = size(y, 1); 
nout = size(Y.P.wav, 2) * nch; 
Y.y = y; 

DCM.Y = Y; % update

%--------------------------------------------------------------------------
% detrend outputs (and inputs)  
% Y.y = spm_detrend(Y.y); 
if DCM.options.centre
        U.u = spm_detrend(U.u); 
end

% generate regressors for confounding effects 
%--------------------------------------------------------------------------
X0 = [];
switch Y.C.type 
    case 'DCT' 
        nf = size(Y.C.period, 1); 
        for i = 1:nf
            period = Y.C.period(i,:); 
            Lorder = fix(2*(v./Y.P.fs)/period(1,2) + 1);
            Horder = fix(2*(v./Y.P.fs)/period(1,1) + 1);
            X0_t = spm_dctmtx(v, Horder);
            X0_t = X0_t(:, Lorder:end);
            X0 = [X0 X0_t];
        end
    case 'User' 
        X0 = [X0 Y.C.X0]; 
end 
X0 = [X0 ones(v, 1)]; 
Y.X0 = X0; 

% calculate sensitivity matrix using precalculated Green's functions 
%--------------------------------------------------------------------------
A = spm_fnirs_sensitivity(DCM); 

%--------------------------------------------------------------------------
% create priors 

%--------------------------------------------------------------------------
% check DCM.d (for nonlinear DCMs)
try
    DCM.options.nonlinear = logical(size(DCM.d,3));
catch
    DCM.d = zeros(n,n,0);
    DCM.options.nonlinear = 0;
end

%--------------------------------------------------------------------------
% specify parameters for spm_int_D (ensuring updates every second or so)
if DCM.options.nonlinear
    M.IS     = 'spm_int_D';
    M.nsteps = round(max(Y.dt,1));
    M.states = 1:n;
else
    M.IS     = 'spm_int';
end

%--------------------------------------------------------------------------
% check for endogenous DCMs, with no exogenous driving effects
if isempty(DCM.c) || isempty(U.u)
    DCM.c  = zeros(n,1);
    DCM.b  = zeros(n,n,1);
    U.u    = zeros(v,1);
    U.name = {'null'};
end
if ~any(spm_vec(U.u)) || ~any(spm_vec(DCM.c))
    DCM.options.stochastic = 1;
end

%--------------------------------------------------------------------------
% priors (and initial states)
[pE,pC,x] = spm_dcm_fnirs_priors(DCM); 

try, M.P     = DCM.options.P;  end      % initial parameters
try, pE      = DCM.options.pE; end      % prior expectation
try, pC      = DCM.options.pC; end      % prior covariance

%--------------------------------------------------------------------------
% eigenvector constraints on pC for large models
if n > DCM.options.nmax % This routine will be skipped. 
    
    %----------------------------------------------------------------------
    % remove confounds and find principal (nmax) modes
    y       = Y.y - Y.X0*(pinv(Y.X0)*Y.y);
    V       = spm_svd(y');
    V       = V(:,1:DCM.options.nmax);
    
    %----------------------------------------------------------------------
    % remove minor modes from priors on A
    j       = 1:(n*n);
    V       = kron(V*V',V*V');
    pC(j,j) = V*pC(j,j)*V';
end

%--------------------------------------------------------------------------
% hyperpriors over precision - expectation and covariance
hE      = sparse(nout,1) + 6;
hC      = speye(nout,nout)/128;
i       = DCM.options.hidden;
hE(i)   = -4;
hC(i,i) = exp(-16);

Y.Q        = spm_Ce(ones(1,nout)*v); % error precision components 

%--------------------------------------------------------------------------
% complete model specification
M.rs = DCM.options.rs;  % radius of source
M.acoef= Y.P.acoef;  % molar absorption coefficients [1/(mM*cm)] 
M.f  = 'spm_fx_fnirs'; % hemodynamic state equation
M.g  = 'spm_gx_state_fnirs'; % optics equation

M.x = x; 
M.A = A; % sensitivity matrix 
M.pE = pE; % prior expectation (parameters)
M.pC = pC; % prior covariance  (parameters)
M.m  = size(U.u,2);
M.n  = size(M.x(:),1);
M.nwav = size(Y.P.wav, 2); 
M.nch = nch; % number of channels

M.l  = nout;
M.hE = hE; % prior expectation (precisions)
M.hC = hC; % prior covariance  (precisions)
M.N  = 32;
M.dt = 16/M.N;
M.ns = v;

%--------------------------------------------------------------------------
% nonlinear system identification (nlsi)
if ~DCM.options.stochastic
    
    %----------------------------------------------------------------------
    % nonlinear system identification (Variational EM) - deterministic DCM
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);
    
    %----------------------------------------------------------------------
    % predicted responses (y) and residuals (R)
    y      = feval(M.IS,Ep,M,U);
    R      = Y.y - y;
    R      = R - Y.X0*spm_inv(Y.X0'*Y.X0)*(Y.X0'*R);
    Ce     = exp(-Eh);   
end

%--------------------------------------------------------------------------
% Bilinear representation and first-order hemodynamic kernel
[M0,M1,L1,L2] = spm_bireduce(M,Ep);
[H0,H1] = spm_kernels(M0,M1,L1,L2,M.N,M.dt);

%--------------------------------------------------------------------------
% and neuronal kernels
L       = sparse(1:n,(1:n) + 1,1,n,length(M0));
[K0,K1] = spm_kernels(M0,M1,L,M.N,M.dt);


%--------------------------------------------------------------------------
% Bayesian inference and variance {threshold: prior mean plus T = 0}
T       = full(spm_vec(pE));
sw      = warning('off','SPM:negativeVariance');
Pp      = spm_unvec(1 - spm_Ncdf(T,abs(spm_vec(Ep)),diag(Cp)),Ep);
Vp      = spm_unvec(full(diag(Cp)),Ep);
warning(sw);
try,  M = rmfield(M,'nograph'); end

%--------------------------------------------------------------------------
% Store parameter estimates
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

%--------------------------------------------------------------------------
% Data ID and log-evidence
if isfield(M,'FS')
    try
        ID = spm_data_id(feval(M.FS,Y.y,M));
    catch
        ID = spm_data_id(feval(M.FS,Y.y));
    end
else
    ID     = spm_data_id(Y.y);
end

%--------------------------------------------------------------------------
% Save approximations to model evidence: negative free energy, AIC, BIC
evidence   = spm_dcm_evidence(DCM);
DCM.F      = F;
DCM.ID     = ID;
DCM.AIC    = evidence.aic_overall;
DCM.BIC    = evidence.bic_overall;

%--------------------------------------------------------------------------
%-Save DCM
if ~isstruct(P)
    save(P, 'DCM', 'F', 'Ep', 'Cp', spm_get_defaults('mat.format')); 
end

if ~nargin
    spm('Pointer','Arrow');
    spm('FigName','Done');
end



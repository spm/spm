function DCM = spm_dcm_tfm(DCM)
% Estimate parameters of a DCM of (induced) cross-spectral density
% FORMAT DCM = spm_dcm_tfm(DCM)
%
% DCM
%    name: name string
%       xY: data   [1x1 struct]
%       xU: design [1x1 struct]
%
%   Sname: cell of source name strings
%       A: {[nr x nr double], [nr x nr double], ...}
%       B: {[nr x nr double], ...}   Connection constraints
%       C: [nr x 1 double]
%
%   options.Nmodes       - number of spatial modes
%   options.Tdcm         - [start end] time window in ms
%   options.Fdcm         - [start end] Frequency window in Hz
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.spatial      - 'ECD', 'LFP' or 'IMG'     (see spm_erp_L)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm.m 4807 2012-07-26 16:15:49Z guillaume $
 
 
% check options
%==========================================================================
drawnow
clear spm_erp_L
name = sprintf('DCM_%s',date);
 
% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                    catch, DCM.name = name;      end
try, Nm    = DCM.options.Nmodes;  catch, Nm = 8;               end
try, onset = DCM.options.onset;   catch, onset     = 60;       end
try, dur   = DCM.options.dur;     catch, dur       = 16;       end
 
% Spatial model
%==========================================================================
model                = 'CMM';
DCM.options.analysis = 'TFR';
DCM.options.Nmodes   = Nm;
DCM.options.model    = model;
DCM.M.dipfit.model   = model;

DCM  = spm_dcm_erp_dipfit(DCM, 1);                  % spatial model
DCM  = spm_dcm_erp_data(DCM);                       % data
Ns   = length(DCM.A{1});                            % number of sources
 

% Design model and exogenous inputs
%==========================================================================
if isempty(DCM.xU.X), DCM.xU.X = sparse(1 ,0); end
if ~isfield(DCM,'C'), DCM.C    = sparse(Ns,0); end
if isempty(DCM.xU.X), DCM.C    = sparse(Ns,0); end

% Neural mass model
%==========================================================================
 
% prior moments on parameters
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
[pE,pC]  = spm_ssr_priors(pE,pC);
 
% initial states and equations of motion
%--------------------------------------------------------------------------
[x,f]    = spm_dcm_x_neural(pE,model);

 
% orders and model
%==========================================================================
nx       = length(spm_vec(x));
nu       = size(pE.C,2);

% create DCM
%--------------------------------------------------------------------------
% DCM.M.FS = 'spm_fs_csd';
DCM.M.IS = 'spm_csd_int';
DCM.M.g  = 'spm_gx_erp';
DCM.M.f  = f;
DCM.M.x  = x;
DCM.M.n  = nx;
DCM.M.m  = nu;
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.hE = 8;
DCM.M.hC = exp(-4);

% solve for steady state
%--------------------------------------------------------------------------
DCM.M.x  = spm_dcm_neural_x(pE,DCM.M);


%-Feature selection using principal components (U) of lead-field
%==========================================================================
DCM.M.U  = spm_dcm_eeg_channelmodes(DCM.M.dipfit,Nm);


% get data-features (in reduced eigenspace)
%--------------------------------------------------------------------------
DCM      = spm_dcm_tfm_data(DCM);

% within-trial effects: adjust onset relative to PST
%--------------------------------------------------------------------------
DCM.M.ons = onset - DCM.xY.pst(1);
DCM.M.dur = dur;
DCM.xU.dt = DCM.xY.dt;

 
 
% complete model specification and invert
%==========================================================================
Nm        = size(DCM.M.U,2);                    % number of spatial modes
Nf        = length(DCM.xY.Hz);                  % number of frequency bins
Nb        = length(DCM.xY.pst);                 % number of time bins
Nt        = length(DCM.xY.csd);                 % number of trial types
DCM.M.l   = Nm;
DCM.M.Hz  = DCM.xY.Hz;
DCM.M.Rft = DCM.xY.Rft;


% precision of noise
%--------------------------------------------------------------------------
Qf        = spm_Q(1/4,Nf,1);
Qt        = spm_Q(1/4,Nb,1);
Q         = kron(Qf,Qt);
DCM.xY.Q  = Q;
DCM.xY.X0 = sparse(Nt*Nm*Nm*Nf*Nb,0);


% adjust gain to accommodate scaling differences among models and data
%==========================================================================


% % check neural activity (without sensor noise) and extrinsic coupling
% %--------------------------------------------------------------------------
% y         = feval(DCM.M.IS,pE,DCM.M,DCM.xU);
% scale     = mean(abs(spm_vec(y)));
% DCM.M.U   = DCM.M.U/sqrt(scale);

[csd,x,x,x,x,x,erp] = spm_csd_int(pE,DCM.M,DCM.xU);
xY.erp = erp;
xY.csd = csd;
spm_dcm_tfm_response(xY,DCM.xY.pst,DCM.xY.Hz)


% Variational Laplace: model inversion
%==========================================================================
[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);


% Data ID
%--------------------------------------------------------------------------
try
    try
        ID  = spm_data_id(feval(DCM.M.FS,DCM.xY.y,DCM.M));
    catch
        ID  = spm_data_id(feval(DCM.M.FS,DCM.xY.y));
    end
catch
    ID  = spm_data_id(DCM.xY.y);
end
 
 
% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
Hc  = spm_csd_mtf(Qp,DCM.M,DCM.xU);                      % prediction
Ec  = spm_unvec(spm_vec(DCM.xY.y) - spm_vec(Hc),Hc);     % prediction error
 
 
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
M             = DCM.M;    
M.dipfit.type = 'LFP';
 
qp        = Qp;
qp.L      = ones(1,Ns);             % set virtual electrode gain to unity
qp.b      = qp.b - 32;              % and suppress non-specific and
qp.c      = qp.c - 32;              % specific channel noise

[Hs,Hz,dtf] = spm_csd_mtf(qp,M,DCM.xU);
[ccf,pst]   = spm_csd2ccf(Hs,DCM.M.Hz);
[coh,fsd]   = spm_csd2coh(Hs,DCM.M.Hz);
DCM.dtf  = dtf;
DCM.ccf  = ccf;
DCM.coh  = coh;
DCM.fsd  = fsd;
DCM.pst  = pst;
DCM.Hz   = Hz;

 
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
 
% and save
%--------------------------------------------------------------------------
DCM.options.Nmodes = Nm;
 
save(DCM.name, 'DCM', spm_get_defaults('mat.format'));

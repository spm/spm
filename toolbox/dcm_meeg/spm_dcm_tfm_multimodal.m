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
%   options.h            - order of (DCT) detrending
%   options.Tdcm         - [start end] time window in ms
%   options.Fdcm         - [start end] Frequency window in Hz
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.spatial      - 'ECD', 'LFP' or 'IMG'     (see spm_erp_L)
%
% Returns:
%
% sensor space
%--------------------------------------------------------------------------
% DCM.csd;                  % conditional cross-spectral density
% DCM.tfm;                  % conditional induced responses
% DCM.dtf;                  % conditional directed transfer functions
% DCM.erp;                  % conditional evoked responses
% DCM.Qu;                   % conditional neuronal responses
% DCM.pst;                  % peristimulus times
% DCM.Hz;                   % frequencies
%
% store estimates in DCM
%--------------------------------------------------------------------------
% DCM.Ep;                   % conditional expectation - parameters
% DCM.Cp;                   % conditional covariance  - parameters
% DCM.Pp;                   % conditional probability - parameters
% DCM.Ce;                   % error covariance
% DCM.F;                    % Laplace log evidence
% DCM.ID;                   % data ID
%
% source space
%--------------------------------------------------------------------------
% DCM.CSD;                  % conditional cross-spectral density
% DCM.TFM;                  % conditional induced responses
% DCM.DTF;                  % conditional directed transfer functions
% DCM.ERP;                  % conditional evoked responses
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm_multimodal.m 6237 2014-10-12 10:08:13Z karl $
 
 
% check options and preliminaries
%==========================================================================
drawnow
clear spm_erp_L
try, DCM = rmfield(DCM,'M'); end
name     = sprintf('DCM_%s',date);
DCM.options.analysis  = 'TFM';
 
% Filename and options
%--------------------------------------------------------------------------
try, name  = DCM.name;            catch, DCM.name = name;  end
try, Nm    = DCM.options.Nmodes;  catch, Nm       = 8;     end
try, onset = DCM.options.onset;   catch, onset    = 60;    end
try, dur   = DCM.options.dur;     catch, dur      = 16;    end
try, h     = DCM.options.h;       catch, h        = 2;     end

 
% Design model and exogenous inputs
%==========================================================================
if isempty(DCM.xU.X), DCM.xU.X = sparse(1,0); end
 
 
% Spatial model
%==========================================================================
model              = 'TFM';
DCM.options.Nmodes = Nm;
DCM.options.model  = model;
DCM.options.Nmax   = 32;
DCM.options.h      = h;
DCM.M.dipfit.model = model;
DCM.M.N            = 2^16;

 
% Get posterior from event-related responses
%==========================================================================
ERP            = DCM;
[pth name]     = fileparts(DCM.name);
ERP.name       = fullfile(pth,[name '_erp']);
ERP.M.TFM      = 1;
ERP.M.hE       = 8;
ERP.M.hC       = 1/128;
ERP            = spm_dcm_erp(ERP);

%-Feature selection using principal components (U) of lead-field
%==========================================================================
try, DCM.M = rmfield(DCM.M,'TFM'); end
clear functions

DCM.xY         = ERP.xY;
DCM.M.dipfit   = ERP.M.dipfit ;
DCM.M.U        = ERP.M.U;
DCM            = spm_dcm_tfm_data(DCM);

% concatenate erp and csd into cell array
%--------------------------------------------------------------------------
DCM.xY.y = {};
for i = 1:numel(DCM.xY.erp)
    DCM.xY.y{end + 1} = DCM.xY.erp{i};
    DCM.xY.y{end + 1} = DCM.xY.csd{i};
end
 
% Use posterior as the prior in a model of induced responses
%==========================================================================

% remove very precise modes from neuronal priors
%--------------------------------------------------------------------------
[u s]     = spm_svd(ERP.Cp);
i         = find(diag(s) >= 1/32);
Cp        = u(:,i)*s(i,i)*u(:,i)';
 
% prior moments on parameters (neuronal and spatial)
%--------------------------------------------------------------------------
pE        = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);
pE        = spm_L_priors(DCM.M.dipfit,pE);
pE        = spm_unvec(spm_vec(ERP.Ep,ERP.Eg),pE);
 
% prior moments on parameters (spectral)
%--------------------------------------------------------------------------
[pE,pC]   = spm_ssr_priors(pE);
pC        = spm_cat(spm_diag({Cp,ERP.Cg,diag(spm_vec(pC))}));
 
% initial states and equations of motion
%--------------------------------------------------------------------------
[x,f,h]   = spm_dcm_x_neural(pE,model);
 
% orders and model
%==========================================================================
nx        = length(spm_vec(x));
nu        = size(pE.C,2);

% hyperpriors (slightly less precise than the ERP)
%--------------------------------------------------------------------------
hE        = 6;
hC        = 1/128;
 
% create DCM
%--------------------------------------------------------------------------
DCM.M.IS  = 'spm_csd_int_IS';
DCM.M.g   = 'spm_gx_erp';
DCM.M.f   = f;
DCM.M.h   = h;
DCM.M.x   = x;
DCM.M.n   = nx;
DCM.M.m   = nu;
DCM.M.pE  = pE;
DCM.M.pC  = pC;
DCM.M.hE  = hE;
DCM.M.hC  = hC;
DCM.M.ns  = length(DCM.xY.pst);
 
% solve for steady state
%--------------------------------------------------------------------------
DCM.M.x   = spm_dcm_neural_x(pE,DCM.M);
 
% within-trial effects: adjust onset relative to PST
%--------------------------------------------------------------------------
DCM.M.ons = onset - DCM.xY.pst(1);
DCM.M.dur = dur;
DCM.xU.dt = DCM.xY.dt;
 
 
% complete model specification and invert
%==========================================================================
Ns        = length(DCM.A{1});                   % number of sources
Nm        = size(DCM.M.U,2);                    % number of spatial modes
Nf        = length(DCM.xY.Hz);                  % number of frequency bins
Nb        = length(DCM.xY.pst);                 % number of time bins
DCM.M.l   = Nm;
DCM.M.Hz  = DCM.xY.Hz;
DCM.M.Rft = DCM.xY.Rft;
DCM.M.pst = DCM.xY.pst;


% precision of noise
%--------------------------------------------------------------------------
Qt        = spm_Q(1/2,Nb,1);
Qf        = spm_Q(1/2,Nf,1);
Qerp      = kron(speye(Nm),Qt);
Qcsd      = kron(speye(Nm*Nm),kron(Qf,Qt));
DCM.xY.Q  = {blkdiag(Qerp,Qcsd*0),blkdiag(Qerp*0,Qcsd)};

Nh        = DCM.options.h;
X0        = spm_dctmtx(Nb,Nh);
DCM.xY.X0 = [kron(speye(Nm),X0); sparse(Nm*Nm*Nf*Nb,Nm*Nh)];


 
% Inspect stability
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if 0
    [csd,erp] = spm_csd_int(pE,DCM.M,DCM.xU);
    xY.erp    = erp;
    xY.csd    = csd;
    spm_figure('GetWin','predicted (a priori)')
    spm_dcm_tfm_response(    xY,DCM.xY.pst,DCM.xY.Hz)
    spm_figure('GetWin','empirical time-frequency responses')
    spm_dcm_tfm_response(DCM.xY,DCM.xY.pst,DCM.xY.Hz)
    
    keyboard
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
% Variational Laplace: model inversion
%==========================================================================
DCM.M.Nmax   = 32;
[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
 
 
% Data ID
%--------------------------------------------------------------------------
ID  = spm_data_id(DCM.xY.y);
 
 
% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%--------------------------------------------------------------------------
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd and erp) - sensor space
%==========================================================================
[csd,erp,tfm,dtf,Hz,pst,x] = spm_csd_int(Qp,DCM.M,DCM.xU);
 
% sensor space
%--------------------------------------------------------------------------
DCM.csd = csd;                  % conditional cross-spectral density
DCM.tfm = tfm;                  % conditional induced responses
DCM.dtf = dtf;                  % conditional directed transfer functions
DCM.erp = erp;                  % conditional evoked responses
DCM.Qu  = x;                    % conditional neuronal responses
DCM.pst = pst;                  % peristimulus times (seconds)
DCM.Hz  = Hz;                   % frequencies
 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM.Ep  = Qp;                   % conditional expectation - parameters
DCM.Cp  = Cp;                   % conditional covariance  - parameters
DCM.Pp  = Pp;                   % conditional probability - parameters
DCM.Ce  = exp(-Eh);             % error covariance
DCM.F   = F;                    % Laplace log evidence
DCM.ID  = ID;                   % data ID
 
% predictions (source space) - cf, a LFP from virtual electrode
%==========================================================================
M             = rmfield(DCM.M,'U');
M.dipfit.type = 'LFP';
clear functions
 
Qp.L    = ones(1,Ns);           % set virtual electrode gain to unity
Qp.b    = Qp.b - 32;            % and suppress non-specific and
Qp.c    = Qp.c - 32;            % specific channel noise
M.pE    = Qp;
M.l     = Ns;

[csd,erp,tfm,dtf] = spm_csd_int(Qp,M,DCM.xU);
 
% source space
%--------------------------------------------------------------------------
DCM.CSD = csd;                  % conditional cross-spectral density
DCM.TFM = tfm;                  % conditional induced responses
DCM.DTF = dtf;                  % conditional directed transfer functions
DCM.ERP = erp;                  % conditional evoked responses

 
save(DCM.name, 'DCM', spm_get_defaults('mat.format'));

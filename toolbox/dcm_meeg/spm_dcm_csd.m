function DCM = spm_dcm_csd(DCM)
% Estimate parameters of a DCM of (complex) cross-spectral density
% FORMAT DCM = spm_dcm_csd(DCM)
%
% DCM
%    name: name string
%       xY: data   [1x1 struct]
%       xU: design [1x1 struct]
%
%   Sname: cell of source name strings
%       A: {[nr x nr double]  [nr x nr double]  [nr x nr double]}
%       B: {[nr x nr double], ...}   Connection constraints
%       C: [nr x 1 double]
%
%   options.Nmodes       - number of spatial modes
%   options.Tdcm         - [start end] time window in ms
%   options.Fdcm         - [start end] Frequency window in Hz
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.spatial      - 'ECD', 'LFP' or 'IMG'     (see spm_erp_L)
%   options.model        - 'ERP', 'SEP', 'CMC', 'LFP', 'NMM' or 'MFM'
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_csd.m 4348 2011-06-10 20:50:23Z karl $
 
 
% check options
%==========================================================================
drawnow
clear spm_erp_L
name = sprintf('DCM_%s',date);
 
% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                      catch, DCM.name = name;      end
try, DCM.name;                      catch, DCM.name = 'DCM_SSR'; end
try, model   = DCM.options.model;   catch, model    = 'NMM';     end
try, spatial = DCM.options.spatial; catch, spatial  = 'LFP';     end
try, Nm      = DCM.options.Nmodes;  catch, Nm = 8;               end
 
% Spatial model
%==========================================================================
DCM.options.Nmodes = Nm;
 
DCM  = spm_dcm_erp_dipfit(DCM, 1);
Ns   = size(DCM.C,1);                                   % number of sources
Nc   = DCM.M.dipfit.Nc;
DCM  = spm_dcm_erp_data(DCM);
DCM.M.dipfit.model = model;
 

% Design model
%==========================================================================
if isempty(DCM.xU.X), DCM.xU.X = sparse(1,0); end

% Neural mass model
%==========================================================================
 
% prior moments on parameters
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);
 
% reduce prior conduction delays for LFP spatial models
%--------------------------------------------------------------------------
if strcmpi(spatial,'LFP')
    DCM.M.Pf.D = [2 16];
end
 
% check to see if neuronal priors have already been specified
%--------------------------------------------------------------------------
try
    if length(spm_vec(DCM.M.pE)) == length(spm_vec(pE));
        pE = DCM.M.pE;
        fprintf('Using existing priors\n')
    end
end
 
% augment with priors on spatial model
%--------------------------------------------------------------------------
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
 
% augment with priors on endogenous inputs (neuronal) and noise
%--------------------------------------------------------------------------
[pE,pC]  = spm_ssr_priors(pE,pC);
 
% initial states and equations of motion
%--------------------------------------------------------------------------
[x,f]    = spm_dcm_x_neural(pE,model);
 
% create DCM
%--------------------------------------------------------------------------
DCM.M.IS = 'spm_csd_mtf';
DCM.M.FS = 'spm_fs_csd';
DCM.M.g  = 'spm_gx_erp';
DCM.M.f  = f;
DCM.M.x  = x;
DCM.M.n  = length(spm_vec(x));
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.hE = 8;
DCM.M.hC = exp(-8);
DCM.M.m  = Ns;
DCM.M.u  = sparse(Ns,1);
 
%-Feature selection using principal components (U) of lead-field
%==========================================================================
 
% Spatial modes
%--------------------------------------------------------------------------
if Nc < Nm
    U     = speye(Nc,Nc);
    DCM.M.U = U;
else
    dGdg  = spm_diff('spm_lx_erp',pE,DCM.M,1);
    L     = spm_cat(dGdg);
    U     = spm_svd(L*L',exp(-8));
    try
        U = U(:,1:Nm);
    end
    DCM.M.U = U;
end
 
% get data-features (in reduced eigenspace)
%--------------------------------------------------------------------------
DCM       = spm_dcm_csd_data(DCM);
 
 
% complete model specification and invert
%==========================================================================
Nm        = size(DCM.M.U,2);                    % number of spatial modes
Nf        = size(DCM.xY.y{1},1);                % number of frequency bins
DCM.M.l   = Nm;
DCM.M.Hz  = DCM.xY.Hz;
 
% precision of noise: AR(1/2)
%--------------------------------------------------------------------------
q         = spm_Q(1/2,Nf,1)*diag(DCM.M.Hz)*spm_Q(1/2,Nf,1);
d         = speye(Nm,Nm);                       % autospectra
Q{1}      = kron(diag(d(:)),q);
d         = 1 - speye(Nm,Nm);                   % crossspectra
Q{2}      = kron(diag(d(:)),q);
DCM.xY.Q  = Q;
DCM.xY.X0 = sparse(Nf,0);


% adjust gain to accommodate scaling differences among models and data
%==========================================================================

% cross-spectral data
%--------------------------------------------------------------------------
y         = spm_vec(DCM.xY.y);
scale     = mean(abs(y));
DCM.xY.y  = spm_unvec(y/scale,DCM.xY.y);

% check neural activity (without sensor noise) and extrinsic coupling
%--------------------------------------------------------------------------
pE.b      = pE.b - 32;
pE.c      = pE.c - 32;
scale     = mean(abs(spm_vec(feval(DCM.M.IS,pE,DCM.M,DCM.xU))));
while scale > 32
    pE.A  = spm_unvec(spm_vec(pE.A) - 1/8,pE.A);
    scale = mean(abs(spm_vec(feval(DCM.M.IS,pE,DCM.M,DCM.xU))));
end
DCM.M.pE.A = pE.A;
DCM.M.U    = DCM.M.U/sqrt(scale)/2;


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

[Hs Hz dtf] = spm_csd_mtf(qp,M,DCM.xU);
[ccf pst]   = spm_csd2ccf(Hs,DCM.M.Hz);
[coh fsd]   = spm_csd2coh(Hs,DCM.M.Hz);
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
 
if spm_check_version('matlab','7') >= 0
    save(DCM.name, '-V6', 'DCM');
else
    save(DCM.name, 'DCM');
end


function DCM = spm_dcm_nfm(DCM)
% Estimate parameters of a DCM of spectral neural field activity
% FORMAT DCM = spm_dcm_nfm(DCM)
%
% DCM
%    name: name string
%       M:  Forward model
%              M.dipfit - lead-field specification
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
%   options.type         - 'ECD', 'LFP' or 'IMG'     (see spm_erp_L)
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_nfm.m 4492 2011-09-16 12:11:09Z guillaume $
 
 
% check options
%==========================================================================
drawnow
clear spm_erp_L
name = sprintf('DCM_%s',date);
 
% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                   catch, DCM.name = name;      end
try, DCM.name;                   catch, DCM.name = 'DCM_SSR'; end
try, model = DCM.options.model;  catch, model    = 'NFM';     end
try, Nm    = DCM.options.Nmodes; catch, Nm       = 8;         end
 
% Spatial model
%==========================================================================
DCM  = spm_dcm_erp_dipfit(DCM, 1);
Ns   = size(DCM.C,1);                                   % number of sources
Nc   = DCM.M.dipfit.Nc;
DCM  = spm_dcm_erp_data(DCM);
DCM.M.dipfit.model = model;
 
 
% Neural mass model
%==========================================================================
 
% prior moments on parameters
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);
 
% augment with priors on spatial model
%--------------------------------------------------------------------------
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
 
% augment with priors on endogenous inputs (neuronal) and noise
%--------------------------------------------------------------------------
[pE,pC]  = spm_ssr_priors(pE,pC);
 
% create DCM
%--------------------------------------------------------------------------
DCM.M.IS = 'spm_nfm_mtf';
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.m  = Ns;
DCM.M.hE = 8;
DCM.M.hC = exp(-8);
 
 
% Spatial modes and data-features
%==========================================================================
DCM       = spm_dcm_ssr_data(DCM);
 
% complete model specification and invert
%==========================================================================
Nt        = size(DCM.xY.y,1);              % number of trials
Nf        = size(DCM.xY.y{1},1);           % number of frequency bins
DCM.M.l   = Nc;
DCM.M.Hz  = DCM.xY.Hz;
 
% precision of noise
%--------------------------------------------------------------------------
DCM.xY.Q  = {spm_Q(1/2,Nf,1)*diag(sqrt(DCM.xY.Hz))*spm_Q(1/2,Nf,1)};
DCM.xY.X0 = sparse(Nf,0);
 
% scale cross-spectral data
%--------------------------------------------------------------------------
y         = spm_vec(DCM.xY.y);
scale     = mean(abs(y));
DCM.xY.y  = spm_unvec(y/scale,DCM.xY.y);
 
% scale neural activity (without sensor noise)
%--------------------------------------------------------------------------
DCM.M.U   = 1;
pE.b      = pE.b - 32;
pE.c      = pE.c - 32;
scale     = mean(abs(spm_vec(feval(DCM.M.IS,pE,DCM.M,DCM.xU))));
DCM.M.U   = DCM.M.U/scale;
 
 
% Model inversion
%--------------------------------------------------------------------------
[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
Ce           = exp(-Eh);
 
% Data ID
%==========================================================================
if isfield(DCM.M,'FS')
    try
        ID  = spm_data_id(feval(DCM.M.FS,DCM.xY.y,DCM.M));
    catch
        ID  = spm_data_id(feval(DCM.M.FS,DCM.xY.y));
    end
else
    ID  = spm_data_id(DCM.xY.y);
end
 
 
% Bayesian inference {threshold = prior} NB Prior on A,B  and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions and error (source space)
%--------------------------------------------------------------------------
Hc  = feval(DCM.M.IS,Qp,DCM.M,DCM.xU);                   % prediction
Ec  = spm_unvec(spm_vec(DCM.xY.y) - spm_vec(Hc),Hc);     % prediction error
 
 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM.Ep = Qp;                   % conditional expectation
DCM.Cp = Cp;                   % conditional covariance
DCM.Pp = Pp;                   % conditional probability
DCM.Hc = Hc;                   % conditional responses (y), channel space
DCM.Rc = Ec;                   % conditional residuals (y), channel space
DCM.Ce = Ce;                   % conditional error covariance
DCM.F  = F;                    % Laplace log evidence
DCM.ID = ID;                   % data ID
 
% and save
%--------------------------------------------------------------------------
DCM.options.Nmodes = Nc;
 
save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
assignin('base','DCM',DCM)

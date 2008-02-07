function DCM = spm_dcm_ssr(DCM)   
% Estimate parameters of a DCM model of cross spectral density of responses
% FORMAT DCM = spm_dcm_ssr(DCM)   
%
% DCM     
%    name: name string
%       M:  Forward model
%              M.dipfit – lead-field specification
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
%   options.type         - 1 - 'ECD (EEG)'
%                          2 - 'ECD (MEG)'
%                          3 - 'Imaging'
%                          4 - 'LFP' 
%                          (see spm_erp_L)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ssr.m 1143 2008-02-07 19:33:33Z spm $
 
 
% check options 
%==========================================================================
clear spm_erp_L
 
% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                  catch, DCM.name           = 'DCM_SSR'; end
try, Nm = DCM.options.Nmodes;   catch, DCM.options.Nmodes = 8; Nm = 8; end
try, h  = DCM.options.h;        catch, h                  = 1;         end
 
% Spatial model
%==========================================================================
DCM     = spm_dcm_erp_dipfit(DCM);
 
% Neural mass model
%==========================================================================
H       = sparse(9,1,1,13,1);   % mixture of states subtending LFP(EEG)
L       = DCM.M.dipfit.L;       % mixture of sources subtending LFP(EEG)
 
% get priors (ensuring C is a diagonal matrix - an input for each source
%--------------------------------------------------------------------------
[pE,pC] = spm_lfp_priors(DCM.A,DCM.B,diag(sum(DCM.C,2)),L,H);
 
 
% create DCM
%--------------------------------------------------------------------------
ns        = size(DCM.C,1);
DCM.M.IS  = 'spm_lfp_mtf';
DCM.M.FS  = 'spm_lfp_sqrt';
DCM.M.f   = 'spm_fx_lfp';
DCM.M.g   = 'spm_gx_lfp';
DCM.M.x   = sparse(ns,13);
DCM.M.n   = ns*13;
DCM.M.pE  = pE;
DCM.M.pC  = pC;
DCM.M.m   = ns;
 
%-Feature selection using principal components (U) of lead-feild
%==========================================================================
 
% Get spatial modes from lead-field
%--------------------------------------------------------------------------
L     = 0;
dGdg  = spm_diff('spm_erp_L',pE,DCM.M,1);
for i = 1:length(dGdg)
    L = L + pC(i,i)*dGdg{i}*dGdg{i}';
end

% eigenspace
%--------------------------------------------------------------------------
U     = spm_svd(L,exp(-8));
try
   DCM.M.U = U(:,1:Nm);
end
Nm    = size(DCM.M.U,2);
 
% get data (in reduced eigen-space) & assume i.i.d noise over CSDs
%--------------------------------------------------------------------------
DCM   = spm_dcm_ssr_data(DCM);


% complete model specification and invert
%========================================================================== 
Nt        = size(DCM.xY.y,1);              % number of trials
Nf        = size(DCM.xY.y{1},1);           % number of frequency bins
DCM.xY.X0 = ones(Nf*Nm*Nm,1);              % confounds (mean in log space)
DCM.xY.Q  = spm_Q(1/2,Nf,1);               % precision of noise AR(1/2)
DCM.M.l   = Nm;
DCM.M.Hz  = DCM.xY.Hz;

% EM: inversion
%--------------------------------------------------------------------------
[Qp,Cp,Ce,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);


 
% Bayesian inference {threshold = prior} NB Prior on A,B  and C = exp(0) = 1
%==========================================================================
warning off
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning on
 
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
DCM.Ce = Ce;                   % ReML error covariance
DCM.F  = F;                    % Laplace log evidence
 
 
% and save
%--------------------------------------------------------------------------
if spm_matlab_version_chk('7.1') >= 0
    save(DCM.name, '-V6', 'DCM');
else
    save(DCM.name, 'DCM');
end
assignin('base','DCM',DCM)
return

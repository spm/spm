function DCM = spm_dcm_ssr(DCM)   
% Estimate parameters of a DCM of cross-spectral density
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
%   options.type         - 'ECD' (1) or 'Imaging' (2) (see spm_erp_L)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ssr.m 1277 2008-03-28 18:36:49Z karl $
 
 
% check options 
%==========================================================================
clear spm_erp_L
 
% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                  catch, DCM.name           = 'DCM_SSR'; end
try, Nm = DCM.options.Nmodes;   catch, DCM.options.Nmodes = 8; Nm = 8; end

% Spatial model
%==========================================================================
DCM  = spm_dcm_erp_data(DCM);
DCM  = spm_dcm_erp_dipfit(DCM);

% Neural mass model
%==========================================================================
 
% get priors (ensuring C is a diagonal matrix - an input for each source
%--------------------------------------------------------------------------
[pE,pC]  = spm_lfp_priors(DCM.A,DCM.B,diag(sum(DCM.C,2)),DCM.M.dipfit);

% create DCM
%--------------------------------------------------------------------------
ns       =  size(DCM.C,1);                         % number of sources
DCM.M.IS = 'spm_lfp_mtf';
DCM.M.FS = 'spm_lfp_sqrt';
DCM.M.f  = 'spm_fx_lfp';
DCM.M.g  = 'spm_gx_erp';
DCM.M.x  = sparse(ns,13);
DCM.M.n  = ns*13;
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.m  = ns;
 
%-Feature selection using principal components (U) of lead-feild
%==========================================================================

% get eigen-space of sources
%--------------------------------------------------------------------------
switch DCM.xY.modality

    case{'EEG','MEG'}

        % Get spatial modes from lead-field
        %------------------------------------------------------------------
        L     = 0;
        dGdg  = spm_diff('spm_erp_L',pE,DCM.M,1);
        for i = 1:length(dGdg)
            L = L + pC(i,i)*dGdg{i}*dGdg{i}';
        end

        % eigenspace
        %------------------------------------------------------------------
        U     = spm_svd(L,exp(-8));
        try
            DCM.M.U = U(:,1:Nm);
        end

    otherwise
        
        % assume a sufficiently small number of LFP channels
        %------------------------------------------------------------------
        DCM.M.U = speye(ns,ns);
end
 
% get data-features (in reduced eigen-space)
%--------------------------------------------------------------------------
DCM      = spm_dcm_ssr_data(DCM);


% complete model specification and invert
%========================================================================== 
Nm       = size(DCM.M.U,2);               % number of spatial modes
Nt       = size(DCM.xY.y,1);              % number of trials
Nf       = size(DCM.xY.y{1},1);           % number of frequency bins
DCM.M.l  = Nm;
DCM.M.Hz = DCM.xY.Hz;

% precision of noise: AR(1/2)
%--------------------------------------------------------------------------
DCM.xY.Q  = spm_Q(1/2,Nf,1);
DCM.xY.X0 = sparse(Nf,0);

% EM: inversion
%--------------------------------------------------------------------------
[Qp,Cp,Ce,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);

 
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
DCM.Ce = Ce;                   % ReML error covariance
DCM.F  = F;                    % Laplace log evidence
 

% and save
%--------------------------------------------------------------------------
DCM.options.Nmodes = Nm;

if spm_matlab_version_chk('7.1') >= 0
    save(DCM.name, '-V6', 'DCM');
else
    save(DCM.name, 'DCM');
end
assignin('base','DCM',DCM)
return

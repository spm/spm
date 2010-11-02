function spm_dcm_post_hoc(P)
% Post hoc optimisation of DCMs (under Laplace approximation)
% FORMAT spm_dcm_post_hoc(P)
%
% P         -  character/cell array of DCM filenames
% name      -  name of DCM output file (will be prefixed by 'DCM_opt_')
%
% This routine searches over all reduced models and uses post-hoc model 
% selection to select the best model. Reduced models here mean all 
% permutations of free parameters (coupling parameters with a non-zero 
% prior covariance), where models are defined in terms of their prior 
% covariance.
% 
% When several DCMs are selected, they are checked to ensure the same free 
% parameters have been specified and the log-evidences are pooled in a 
% fixed effects fashion.
% 
% The outputs of this routine are graphics reporting the model reduction 
% (optimisation) and a DCM_opt_??? for every input DCM that contains the 
% reduced conditional parameters estimates (for simplicity, the original 
% kernels and predicted states are retained). The structural and function 
% (spectral embedding) graphs are based on are based Bayesian parameter 
% averages over multiple DCMs.
%  
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_post_hoc.m 4108 2010-11-02 20:24:02Z karl $
 
% get filenames
%--------------------------------------------------------------------------
try
    P;
catch
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
 
if ischar(P), P = cellstr(P); end
N = numel(P);
 
%-Check models are compatible in terms of their prior variances
%==========================================================================
for j = 1:N
    
    % get prior covariances
    %----------------------------------------------------------------------
    load(P{j});
    pC    = diag(DCM.M.pC);
    
    % and compare it with the first model
    %----------------------------------------------------------------------
    if j == 1
        C = pC;
    else
        if any(~pC - ~C)
            fprintf('Please check model %i for compatibility',model)
            return
        end
    end
end
 
% find free coupling parameters
%--------------------------------------------------------------------------
k     = spm_fieldindices(DCM.Ep,'A','B','D');
k     = k(~~C(k));
 
% Create model space in terms of parameter indices
%--------------------------------------------------------------------------
K     = spm_perm_mtx(length(k));
n     = length(C);
 
%-Loop through models and get log-evidences
%==========================================================================
 
for j = 1:N
    
    load(P{j});
    
    % Get priors and posteriors
    % ---------------------------------------------------------------------
    qE    = DCM.Ep;
    qC    = DCM.Cp;
    pE    = DCM.M.pE;
    pC    = DCM.M.pC;
    
    % model search over new prior (covariance)
    % ---------------------------------------------------------------------
    for i = 1:length(K)
        R      = speye(n,n) - sparse(k,k,K(i,:),n,n);
        G(i,j) = spm_log_evidence(qE,qC,pE,pC,pE,R*pC*R);
    end
    
end
 
% Pooled model evidence
% =========================================================================
S     = sum(G,2);
S     = S - S(end);
p     = exp(S - max(S));
p     = p/sum(p);
 
% Get selected model (prior covariance)
%--------------------------------------------------------------------------
[q i] = max(p);
R     = speye(n,n) - sparse(k,k,K(i,:),n,n);
rC    = R*pC*R;
 
 
% Show results
% -------------------------------------------------------------------------
spm_figure('Getwin','Graphics'); clf
 
subplot(2,2,1)
if length(K) > 32, plot(S,'k'), else, bar(S,'c'), end
title('log-posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('log-probability','FontSize',12)
axis square
 
subplot(2,2,2)
if length(K) > 32, plot(p,'k'), else, bar(p,'r'), end
title('model posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('probability','FontSize',12)
axis square
 
 
% conditional estimates of selected model
% =========================================================================
Eq    = 0;
Pq    = 0;
for j = 1:N
    
    % Get priors and posteriors
    % ---------------------------------------------------------------------
    load(P{j});
    
    qE   = DCM.Ep;
    qC   = DCM.Cp;
    pE   = DCM.M.pE;
    pC   = DCM.M.pC;
    
    % Get posterior of selected model - rC
    % ---------------------------------------------------------------------
    [F Ep Cp] = spm_log_evidence(qE,qC,pE,pC,pE,rC);
    
    % Bayesian parameter average
    %----------------------------------------------------------------------
    Pp   = spm_inv(Cp);
    Eq   = Eq + Pp*spm_vec(Ep);
    Pq   = Pq + Pp;
 
 
    % Put reduced conditional estimates in DCM
    % =====================================================================
 
    % Bayesian inference and variance
    %----------------------------------------------------------------------
    sw   = warning('off','SPM:negativeVariance');
    Pp   = spm_unvec(1 - spm_Ncdf(DCM.T,abs(spm_vec(Ep)),diag(Cp)),Ep);
    Vp   = spm_unvec(diag(Cp),Ep);
    warning(sw);
    
    
    % Store parameter estimates
    %----------------------------------------------------------------------
    DCM.M.pC = rC;
    DCM.Ep   = Ep;
    DCM.Cp   = Cp;
    DCM.Pp   = Pp;
    DCM.Vp   = Vp;   
    
    % Save approximations to model evidence: negative free energy, AIC, BIC
    %----------------------------------------------------------------------
    evidence = spm_dcm_evidence(DCM);
    DCM.F    = F;
    DCM.AIC  = evidence.aic_overall;
    DCM.BIC  = evidence.bic_overall;
    
    
    %-Save optimised DCM
    %======================================================================
    [pth, name] = fileparts(P{j});
    filename    = fullfile(pth,['DCM_opt_' name(4:end)]);
    if spm_matlab_version_chk('7') >= 0
        save(filename, 'DCM', '-V6');
    else
        save(filename, 'DCM');
    end
    
end
 
 
 
% Show full and reduced conditional estimates (for last DCM)
%--------------------------------------------------------------------------
spm_figure('Getwin','Graphics');
 
i     = spm_fieldindices(DCM.Ep,'A','B','C','D');
qE    = spm_vec(qE);
Ep    = spm_vec(Ep);
 
subplot(2,2,3)
spm_plot_ci(qE(i),qC(i,i))
title('MAP connections (full)','FontSize',16)
axis square
a   = axis;
 
subplot(2,2,4)
spm_plot_ci(Ep(i),Cp(i,i))
title('MAP connections (reduced)','FontSize',16)
axis square
axis(a)
 
% Show structural and functional graphs
%--------------------------------------------------------------------------
spm_figure('Getwin','Graph'); clf
 
% Bayesian parameter average
%--------------------------------------------------------------------------
Cq     = spm_inv(Pq);
Eq     = Cq*Eq;
Eq     = spm_unvec(Eq,pE);
A      = Eq.A + sum(Eq.B,3) + sum(Eq.D,3);
 
spm_dcm_graph(DCM.xY,A)


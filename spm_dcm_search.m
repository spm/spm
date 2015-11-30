function spm_dcm_search(P)
% Post hoc optimisation of DCMs (under Laplace approximation)
% FORMAT spm_dcm_search(P)
%
% P         -  character/cell array of DCM filenames
%
%--------------------------------------------------------------------------
% spm_dcm_search operates on different DCMs of the same data to identify
% the best model. It will invert the full model whose free-parameters are
% the union (superset) of all free parameters in each model specified. The
% routine then uses a post hoc selection procedure to evaluate the log-
% evidence and conditional density over free-parameters of each model
% specified.
%
% The DCM specified does not need to be estimated. spm_dcm_search will 
% invert the requisite (full DCM) automatically.
%
% The outputs of this routine are graphics reporting the model space search
% (optimisation) and a DCM_optimum (in the first DCMs directory) for the
% best DCM. The structural and function (spectral embedding) graphs are
% based on this DCM.
%
% DCM_optimum  contains the fields:
%        DCM.P   - character/cell array of DCM filenames
%        DCM.PF  - their associated free energies
%        DCM.PP  - and posterior (model) probabilities
%
% In addition, the free energies and posterior estimates of each DCM in P 
% are saved for subsequent searches over different partitions of model 
% space.
%
% See alos: spm_dcm_post_hoc.m
%
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_search.m 6615 2015-11-30 12:56:02Z peter $
 
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
 
%-Check models are compatible in terms of their data
%==========================================================================
for j = 1:N
    
    % get prior covariances
    %----------------------------------------------------------------------
    load(P{j});
    
    % and compare it with the first model
    %----------------------------------------------------------------------
    if j == 1
        Y = DCM.Y.y;
    else
        try
            if any(any(Y - DCM.Y.y))
                fprintf('Please check model %i for compatibility',j)
                return
            end
        catch
            fprintf('Please check model %i for compatibility',j)
            return
        end
    end
    
    % accumate model in terms of which parameters are free
    %----------------------------------------------------------------------
    A = DCM.a;
    B = DCM.b;
    C = DCM.c;
    D = DCM.d;
    
    % Get full models free parameters
    %----------------------------------------------------------------------
    if j == 1
        a    = A;
        b    = B;
        c    = C;
        d    = D;
    else
        a    = a | A;
        b    = b | B;
        c    = c | C;
        d    = d | D;
    end
end
 
%-Estimate full model
%==========================================================================
DCM.a = a;
DCM.b = b;
DCM.c = c;
DCM.d = d;
 
DCM.name = 'optimum';
 
% Get full priors and posteriors
% -------------------------------------------------------------------------
FUL   = spm_dcm_estimate(DCM);
 
qE    = FUL.Ep;
qC    = FUL.Cp;
pE    = FUL.M.pE;
pC    = FUL.M.pC;
 
%-Loop through models and get log-evidences
%==========================================================================
for j = 1:N
    
    load(P{j});
    
    % Fix for endogenous DCM (to match spm_dcm_estimate)
    % ---------------------------------------------------------------------
    if isempty(DCM.c) || isempty(U.u)
        DCM.c  = zeros(DCM.n,1);
        DCM.b  = zeros(DCM.n,DCM.n,1);
    end
    
    % Get model (priors) and evaluate (reduced) free-energy and posteriors
    % ---------------------------------------------------------------------
    [rE,rC]   = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);
    [F,Ep,Cp] = spm_log_evidence(qE,qC,pE,pC,rE,rC);
    
    % Put reduced conditional estimates in DCM
    % =====================================================================
    
    % Bayesian inference and variance
    %----------------------------------------------------------------------
    sw       = warning('off','SPM:negativeVariance');
    Pp       = spm_unvec(1 - spm_Ncdf(0,abs(spm_vec(Ep)),diag(Cp)),Ep);
    Vp       = spm_unvec(diag(Cp),Ep);
    warning(sw);
    
    
    % Store parameter estimates
    %----------------------------------------------------------------------
    DCM.M.pC = rC;
    DCM.Ep   = Ep;
    DCM.Cp   = Cp;
    DCM.Pp   = Pp;
    DCM.Vp   = Vp;
    DCM.T    = 0;
    
    % Store predictions of states from full model for simplicity
    %----------------------------------------------------------------------
    DCM.Ce   = FUL.Ce;
    DCM.H1   = FUL.H1;
    DCM.K1   = FUL.K1;
    DCM.R    = FUL.R;
    DCM.y    = FUL.y;
    
    % Save approximations to model evidence: negative free energy, AIC, BIC
    %----------------------------------------------------------------------
    evidence = spm_dcm_evidence(DCM);
    DCM.F    = F;
    DCM.AIC  = evidence.aic_overall;
    DCM.BIC  = evidence.bic_overall;
    
    % Save DCM
    %======================================================================
    save(P{j},'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));
    
    % Record free-energy
    %----------------------------------------------------------------------
    G(j) = F;
    
end
 
% Model evidences and best model
% =========================================================================
G     = G - min(G);
p     = exp(G - max(G));
p     = p/sum(p);
 
% Get selected model
%--------------------------------------------------------------------------
[q,j] = max(p);
load(P{j});
 
i   = spm_fieldindices(DCM.Ep,'A','B','C','D');
qE  = spm_vec(FUL.Ep);
Ep  = spm_vec(DCM.Ep);
qC  = DCM.Cp;
Cp  = DCM.Cp;
F   = DCM.F;
 
% Show results
% -------------------------------------------------------------------------
spm_figure('Getwin','Graphics'); clf
 
subplot(2,2,1)
if length(P) > 32, plot(G,'k'), else bar(G,'c'), end
title('log-posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('log-probability','FontSize',12)
axis square
 
subplot(2,2,2)
if length(P) > 32, plot(p,'k'), else, bar(p,'r'), end
title('model posterior','FontSize',16)
xlabel('model','FontSize',12)
ylabel('probability','FontSize',12)
axis square
 
% Show full and reduced conditional estimates (for optimum DCM)
%--------------------------------------------------------------------------
subplot(2,2,3)
spm_plot_ci(qE(i),qC(i,i))
title('MAP connections (full)','FontSize',16)
axis square
a   = axis;
 
subplot(2,2,4)
spm_plot_ci(Ep(i),Cp(i,i))
title('MAP connections (optimum)','FontSize',16)
axis square
axis(a)
 
% Show structural and functional graphs
%--------------------------------------------------------------------------
spm_figure('Getwin','Graph'); clf
 
spm_dcm_graph(DCM.xY,DCM.Ep.A)
 
 
%-Save optimum and full DCM
%==========================================================================
DCM.P  = P;
DCM.PF = G;
DCM.PP = p;
 
% Reduced model (optimum)
%--------------------------------------------------------------------------
pth      = fileparts(P{1});
filename = fullfile(pth,'DCM_optimum.mat');
save(filename,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));
 
% Full model
%--------------------------------------------------------------------------
DCM      = FUL;
Ep       = FUL.Ep;
Cp       = FUL.Cp;
F        = FUL.F;
filename = fullfile(pth,'DCM_full');
save(filename,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));

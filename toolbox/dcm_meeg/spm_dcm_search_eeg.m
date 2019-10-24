function [DCM,RCM]   = spm_dcm_search_eeg(P,SAVE_DCM)
% Bayesian model reduction (under Laplace approximation)
% FORMAT [DCM,RCM] = spm_dcm_search_eeg(P,SAVE_DCM)
%
% P - {Nsub x Nmodel} cell array of DCM filenames or model structures for 
%      Nsub subjects, where each model is reduced independently for each
%      subject
%      
% SAVE_DCM - optional flag to save DCMs
% 
% DCM - reduced (best) DCM
% RCM - reduced DCM array
%
% Each reduced model requires DCM.A,DCM.B,DCM.C and DCM.options.model
% or the implicit prior expectation and covariances in DCM.pE and DCM.pC
% if the reduce models are specified explicitly in terms of prior
% expectations and covariances (pE and pC) these will be used first.
%
%--------------------------------------------------------------------------
% spm_dcm_search_eeg operates on different DCMs of the same data to find
% the best model. It assumes the full model - whose free-parameters are
% the union (superset) of all free parameters in each model - has been
% inverted. A post hoc selection procedure is used to evaluate the log-
% evidence and conditional density over free-parameters of each model
% specified.
%
% Reduced models can be specified either in terms of the allowable
% connections (specified in the DCM.A, DCM.B and DCM.C fields) or the
% resulting prior density (specified in DCM.pE and DCM.pC).  If the
% latter exist, they will be used as the model specification.
%
% The outputs of this routine are graphics reporting the model space search
% (optimisation) and a DCM_optimum (in the first DCMs directory) for the
% best DCM. The structural and function (spectral embedding) graphs are
% based on this DCM.
%
% Conditional esimates (Ep, Cp and F values) in DCM_??? (specifed by P) are
% replaced by their reduced estimates (but only these estimates) in rDCM_???
%
% DCM_optimum (saved with nargout = 0) contains the fields:
%
%        DCM.Pname - character/cell array of DCM filenames
%        DCM.PF    - their associated free energies
%        DCM.PP    - and posterior (model) probabilities
%
% If requested, the free energies and posterior estimates of each DCM in P
% are saved for subsequent searches over different partitions of model
% space.
%
% See also: spm_dcm_post_hoc.m, spm_dcm_group and spm_dcm_bma
%
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_search_eeg.m 7679 2019-10-24 15:54:07Z spm $

% get filenames and set up
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if nargin < 2,  SAVE_DCM = 0;    end
if ischar(P),   P = cellstr(P);  end
if isstruct(P), P = {P};         end

% number of subjects and models: BMR over models (rows) for each subject
%--------------------------------------------------------------------------
[Ns N] = size(P);

if Ns > 2
    for i = 1:Ns
        [p,q]    = spm_dcm_search_eeg(P(i,:));
        DCM{i,1} = p;
        RCM(i,:) = q;
    end
    return
end

%-Check models are compatible in terms of their data
%==========================================================================

% establish whether reduced models are specifed in term of pE,pC or A.B.C
%--------------------------------------------------------------------------
try
    
    % number of free parameters
    %----------------------------------------------------------------------
    for j = 1:N
        try, load(P{j}); catch, DCM = P{j}; end
        par(:,j) = logical(diag(DCM.M.pC));
    end
    
    
    
catch
    
    % number of free parameters
    %----------------------------------------------------------------------
    for j = 1:N
        try, load(P{j}); catch, DCM = P{j}; end
        try
            par(:,j) = logical(spm_vec(DCM.A,DCM.B,DCM.C));
        catch
            par(:,j) = logical(spm_vec(DCM.a,DCM.b,DCM.c));
        end
        
    end
end

%-load full model
%==========================================================================
[i j] = max(sum(par));
try, load(P{j}); catch, DCM = P{j}; end

% check this is a full model and that is has been inverted
% -------------------------------------------------------------------------
if any(any(par,2) > par(:,j))
    fprintf('\nPlease ensure your models are nested\n')
    return
end

if ~all(isfield(DCM,{'Ep','Cp'}))
    fprintf('\nPlease invert model %i\n',j)
    return
end


% directory to save reduced [r]DCMs
% -------------------------------------------------------------------------
pathname = spm_file(DCM.name,'path');
if ~exist(pathname,'dir');
    pathname = pwd;
end


% Get full priors and posteriors
% -------------------------------------------------------------------------
options = DCM.options;
qE      = DCM.Ep;
qC      = DCM.Cp;
pE      = DCM.M.pE;
pC      = DCM.M.pC;

%-Loop through models and get log-evidences
%==========================================================================
name  = cell(N,1);
G     = zeros(N,1);
for j = 1:N
    
    % Get reduced model specification
    % ---------------------------------------------------------------------
    try, load(P{j}); catch, DCM = P{j}; end
    
    % Get model (priors)
    % ---------------------------------------------------------------------
    try
        rE   = DCM.M.pE;
        rC   = DCM.M.pC;

    catch
        
        % get priors from alterntie model specification
        %------------------------------------------------------------------
        if isfield(options,'analysis')
            if strcmpi(options.analysis,'IND')
                [rE,gE,rC] = spm_ind_priors(DCM.A,DCM.B,DCM.C,Nf);
            else
                [rE,rC] = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,options.model);
            end
        else
            [rE,rC] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,options);
        end
    end
    
    
    % evaluate (reduced) free-energy and posteriors
    % ---------------------------------------------------------------------
    [F,Ep,Cp] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
    
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
    DCM.M.pE = rE;
    DCM.M.pC = rC;
    DCM.Ep   = Ep;
    DCM.Cp   = Cp;
    DCM.Pp   = Pp;
    DCM.Vp   = Vp;
    DCM.F    = F;
    
    filename = spm_file(DCM.name,'basename');
    name{j}  = filename(5:end);
    filename = spm_file(DCM.name,'filename');
    filename = spm_file(filename,'prefix','r');
    filename = fullfile(pathname,filename);
    
    
    % Save DCM
    %======================================================================
    if SAVE_DCM
        save(filename,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));
    end
    if nargout > 1
        RCM{j} = DCM;
    end
    
    % Record free-energy and MAP estimates
    %----------------------------------------------------------------------
    G(j)     = F;
    
end

% Model evidences and best model
% =========================================================================
G     = G - max(G);
p     = exp(G - max(G));
p     = p/sum(p);
[q,j] = max(p);

% Get selected model
%--------------------------------------------------------------------------
try, load(P{j}); catch, DCM = P{j}; end


qE    = spm_vec(qE);
Ep    = spm_vec(DCM.Ep);
Cp    = DCM.Cp;
F     = DCM.F;
try
    i = find(diag(pC));
catch
    i = find(spm_vec(pC));
end
j     = spm_fieldindices(pE,'A','B','C');
i     = j(ismember(j,i));

% Show results: Graphics
%==========================================================================
spm_figure('Getwin','Graphics'); clf

subplot(2,2,1)
if length(P) > 32, plot(G,'k'), else, bar(diag(G),N), legend(name), end
title('log-posterior','FontSize',16)
xlabel('model','FontSize',12), ylabel('log-probability','FontSize',12)
axis square

subplot(2,2,2)
if length(P) > 32, plot(p,'k'), else, bar(diag(p),N), end
title('model posterior','FontSize',16)
xlabel('model','FontSize',12), ylabel('probability','FontSize',12)
axis square

% Show full and reduced conditional estimates (for optimum DCM)
%--------------------------------------------------------------------------
subplot(2,2,3)
spm_plot_ci(qE(i),qC(i,i))
title('MAP connections (full)','FontSize',16)
axis square, a = axis;

subplot(2,2,4)
spm_plot_ci(Ep(i),Cp(i,i))
title('MAP connections (optimum)','FontSize',16)
axis square, axis(a)



%-Save best DCM
%==========================================================================
DCM.Pname = name;
DCM.PF    = G;
DCM.PP    = p;

% Reduced model (optimum)
%--------------------------------------------------------------------------
if SAVE_DCM
    filename = fullfile(pathname,'DCM_optimum.mat');
    save(filename,'DCM','F','Ep','Cp', spm_get_defaults('mat.format'));
end



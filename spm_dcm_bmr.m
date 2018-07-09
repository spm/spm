function [RCM,BMC,BMA] = spm_dcm_bmr(P,field)
% Bayesian model reduction (under Laplace approximation)
% FORMAT [RCM,BMC,BMA] = spm_dcm_bmr(P,[field])
%
% P     - {Nsub x Nmodel} cell array of DCM filenames or model structures  
%         of Nsub subjects, where each model is reduced independently 
%
% field - parameter fields in DCM{i}.Ep to plot, or the fields to search 
%         if only one DCM is provided per subject [default: {'A','B'}]
%      
% RCM   - reduced DCM array
% BMC   - (Nsub) summary structure 
%          BMC.name - character/cell array of DCM filenames
%          BMC.F    - their associated free energies
%          BMC.P    - and posterior (model) probabilities
% BMA   - Baysian model average (see spm_dcm_bma)
%__________________________________________________________________________
% 
% spm_dcm_bmr operates on different DCMs of the same data (rows) to find
% the best model. It assumes the full model - whose free-parameters are
% the union (superset) of all free parameters in each model - has been
% inverted. A post hoc selection procedure is used to evaluate the log-
% evidence and conditional density over free-parameters of each model
% specified.
%
% Reduced models can be specified either in terms of the allowable
% connections (specified in the DCM.A/a, DCM.B/b and DCM.C/c fields) or the
% resulting prior density (specified in DCM.pE and DCM.pC).  If the
% latter exist, they will be used as the model specification.
%
% If a single subject (DCM) is specified, an exhaustive search will
% be performed.
%
% The outputs of this routine are graphics reporting the model space search
% (optimisation) and the reduced (cell array of) DCM structures.
%
% See also: spm_dcm_post_hoc.m, spm_dcm_bpa, spm_dcm_peb and spm_dcm_bma
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_bmr.m 7369 2018-07-09 09:58:03Z peter $


% get filenames and set up
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P),   P = cellstr(P);  end
if isstruct(P), P = {P};         end

% fields to plot
%--------------------------------------------------------------------------
if nargin < 2;
    field = {'A','B'};
end
if ischar(field)
    field = {field};
end

% number of subjects and models: BMC over models (rows) for each subject
%--------------------------------------------------------------------------
[Ns,N] = size(P);
if Ns > 1
    for i = 1:Ns
        if nargout < 3
            [p,q]    = spm_dcm_bmr(P(i,:),field);
        else
            [p,q,r]  = spm_dcm_bmr(P(i,:),field);
            BMA(i)   = r;            
        end
        RCM(i,:) = p;
        BMC(i)   = q;        
    end
    return
end

% exhaustive search
%--------------------------------------------------------------------------
if N < 2
    [RCM,BMC,BMA] = spm_dcm_bmr_all(P{1},field);
    return
end


%-Check models are compatible in terms of their data
%==========================================================================

% number of free parameters
%----------------------------------------------------------------------
for j = 1:N
    try, load(P{j}); catch, DCM = P{j}; end
    [i,rC,rE,Np] = spm_find_pC(DCM);
    par(:,j)     = sparse(i,1,true,Np,1);
end

%-Load full model
%==========================================================================
[i,j] = max(sum(par));
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

% Get full priors and posteriors
%--------------------------------------------------------------------------
try, DCMname = spm_file(DCM.name,'basename'); catch, DCMname = 'DCM'; end

qE    = DCM.Ep;
qC    = DCM.Cp;
pE    = DCM.M.pE;
pC    = DCM.M.pC;

%-Loop through models and get log-evidences
%==========================================================================
name  = cell(N,1);
RCM   = cell(1,N);
G     = zeros(N,1);
for j = 1:N
    
    % Get reduced model specification
    % ---------------------------------------------------------------------
    try, load(P{j},'DCM'); catch, DCM = P{j}; end
    
    % Get model (priors)
    % ---------------------------------------------------------------------
    [i,rC,rE] = spm_find_pC(DCM);
    
    % evaluate (reduced) free-energy and posteriors
    %----------------------------------------------------------------------
    [F,Ep,Cp] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
    
    % Put reduced conditional estimates in DCM
    %======================================================================
    
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

    % Save DCM and record free-energy
    %----------------------------------------------------------------------
    name{j}  = ['BMC_' DCMname];
    RCM{j}   = DCM;
    G(j)     = DCM.F;
    
end

% Model evidences and best model
%==========================================================================
G     = G - max(G);
p     = exp(G - max(G));
p     = p/sum(p);

%-summary structure
%--------------------------------------------------------------------------
BMC.name = name;
BMC.F    = G;
BMC.P    = p;

% Get and display selected model
%==========================================================================
qE    = spm_vec(qE);
i     = spm_find_pC(pC,pE,field);

if nargout > 2
    
    % Baysian model average
    %----------------------------------------------------------------------
    BMA = spm_dcm_bma(RCM);
    str = 'MAP (average)';
    Ep  = spm_vec(BMA.Ep);
    Cp  = spm_vec(BMA.Cp);
    
else
    
    % Baysian model selection
    %----------------------------------------------------------------------
    [q,j] = max(p);
    str   = 'MAP (selected)';
    Ep    = spm_vec(RCM{j}.Ep);
    Cp    = diag(RCM{j}.Cp);
    
end

if strcmp(field{1},'none'), return, end


% Show results: Graphics
%--------------------------------------------------------------------------
spm_figure('Getwin','Graphics'); clf

subplot(2,2,1)
if length(P) > 32, plot(G,'k'), else, bar(diag(G),N), legend(name), end
title('log-posterior','FontSize',16)
xlabel('model','FontSize',12), ylabel('log-probability','FontSize',12)
axis square

subplot(2,2,2), bar(p,'k')
title('model posterior','FontSize',16)
xlabel('model','FontSize',12), ylabel('probability','FontSize',12)
axis square

% Show full and reduced conditional estimates (for optimum DCM)
%--------------------------------------------------------------------------
subplot(2,2,3), spm_plot_ci(qE(i),qC(i,i))
title('MAP (full)','FontSize',16)
xlabel('parameter'), ylabel('size')
axis square, a = axis;

subplot(2,2,4), spm_plot_ci(Ep(i),Cp(i))
title(str,'FontSize',16)
xlabel('parameter'), ylabel('size')
axis square, axis(a)
drawnow

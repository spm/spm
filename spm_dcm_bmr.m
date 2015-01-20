function [RCM,BMR] = spm_dcm_bmr(P)
% Bayesian model reduction (under Laplace approximation)
% FORMAT [RCM,BMR] = spm_dcm_bmr(P)
%
% P   - {Nsub x Nmodel} cell array of DCM filenames or model structures for 
%       Nsub subjects, where each model is reduced independently for each
%       subject
%      
% RCM - reduced DCM array
% BMR - (Nsub) summary structure 
%        BMR.name - character/cell array of DCM filenames
%        BMR.F    - their associated free energies
%        BMR.P    - and posterior (model) probabilities
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
% The outputs of this routine are graphics reporting the model space search
% (optimisation) and the reduced (cell array of) DCM structures.
%
% See also: spm_dcm_post_hoc.m, spm_dcm_bpa, spm_dcm_peb and spm_dcm_bma
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_bmr.m 6309 2015-01-20 21:01:36Z spm $


% get filenames and set up
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P),   P = cellstr(P);  end
if isstruct(P), P = {P};         end

% number of subjects and models: BMR over models (rows) for each subject
%--------------------------------------------------------------------------
[Ns,N] = size(P);

if Ns > 2
    for i = 1:Ns
        [p,q]    = spm_dcm_bmr(P(i,:));
        RCM(i,:) = p;
        BMR(i)   = q;
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
try, options = DCM.options;                   catch, options = {};    end
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
    name{j}  = ['BMR_' DCMname];
    RCM{j}   = DCM;
    G(j)     = DCM.F;
    
end

% Model evidences and best model
%==========================================================================
G     = G - max(G);
p     = exp(G - max(G));
p     = p/sum(p);

%-Save summary structure
%--------------------------------------------------------------------------
BMR.name = name;
BMR.F    = G;
BMR.P    = p;

% Get and display selected model
%==========================================================================
[q,j] = max(p);
try, load(P{j}); catch, DCM = P{j}; end

qE    = spm_vec(qE);
Ep    = spm_vec(DCM.Ep);
Cp    = DCM.Cp;
try
    i = find(diag(pC));
catch
    i = find(spm_vec(pC));
end
j     = spm_fieldindices(pE,'A','B','C');
i     = j(ismember(j,i));

% Show results: Graphics
%--------------------------------------------------------------------------
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

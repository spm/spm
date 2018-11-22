function [d,BMA,PEBs] = spm_dcm_bdc(GCMs,field,M,ynames,models,noplot)
% Compare datasets using DCM and PEB (Bayesian data comparison)
%
% Performs the following procedure:
%
% 1. Identifies the optimum reduced model (PEB) collapsed over datasets.
% 2. For each dataset, creates reduced GCMs based on the optimum model
%    above and estimates a per-dataset PEB model.
% 3. Computes a series of measures on each dataset's PEB model
%
% Inputs:
%
% GCMs   - {1 x ny} cell array of filenames, each of which is a GCM file.
% field  - {1 x nf} cell array of field names, e.g. {'A','B'}
% M      - (Optional) struct for configuring PEB. See spm_dcm_peb.m .
% M.gamma- (Optional) prior variance for switched off parameters when 
%          automatically forming a model space for KL-over-models.
% M.bmr  - (Optional) if true, searches for an optimal model across all
%          datasets before performing BDC. If false, all parameters from
%          the DCM are included. [default: true];
% ynames - (Optional) {1 x ny} cell array of names for each dataset
% models - (Optional) model space to use for computing the information gain
%          over models. Accepts a nested cell array of parameter names, or
%          a binary matrix (models x parameters) with which parameters
%          to switch on or off in each model.
% noplot - (Optional) if true, does not show plots.
%
% ny = number of datasets, nf = number of DCM fields, nm = number of models
% 
% Returns:
% 
% d.precisions - [np x ny] Precision of each DCM parameter in each dataset
% d.dcm_negent - [1 x ny]  Negative entropy (certainty) of DCM parameters
% d.rfx_negent - [1 x ny]  Negative entropy (certainty) of the estimated
%                          between-subject variability
% d.complexity - [1 x ny]  Number of effective parameters in the model
% d.model_F    - [nm x ny] Free energy of each candidate model
% d.model_P    - [nm x ny] Posterior probability of each candidate model
% d.model_KL   - [1 x ny]  Ability to disriminate between similar models
%
% BMA        - Bayesian model average across all datasets
% PEBs       - [1 x ny] PEB for each dataset
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
 
% Peter Zeidman & Samira Kazan
% $Id: spm_dcm_bdc.m 7495 2018-11-22 15:52:07Z peter $
 
if nargin < 5
    models = [];
end

if nargin < 6
    noplot = false;
end

% Load and concatenate models across datasets
% -------------------------------------------------------------------------
GCM = {};
for i = 1:length(GCMs)    
    if ischar(GCMs{i})
        gcm = load(GCMs{i});
        gcm = gcm.GCM;
    else
        gcm = GCMs{i};
    end
    
    GCM = [GCM; gcm(:,1)];
end

ns = size(gcm,1);   % Number of subjects
ny = length(GCMs);  % Number of datasets
nm = ns*ny;         % Number of DCMs

% Set defaults
% -------------------------------------------------------------------------
if nargin < 2 || isempty(field), field = {'B'}; end
if nargin < 3, M = []; end
if nargin < 4 || isempty(ynames)
    ynames = cell(1,ny); 
    for i = 1:ny
        ynames{i} = sprintf('Dataset %d',i); 
    end
end

if isempty(M)
    M       = struct();
    M.alpha = 1;
    M.beta  = 16;
    M.hE    = 0;
    M.hC    = 1/16;
    M.Q     = 'single';
    M.X     = ones(nm,1);
end

try
    bmr = M.bmr;
catch
    bmr = true;
end

try
    gamma = M.gamma;
catch
    gamma = 0; % Point null hypothesis
end

% Identify optimum model structure across all datasets
% -------------------------------------------------------------------------
if bmr
    % Perform PEB analysis
    PEB = spm_dcm_peb(GCM,M,field);

    % Prune full model
    BMA = spm_dcm_peb_bmc(PEB);

    % Identify disabled PEB parameters
    is_disabled = find(BMA.Pp < 1/1000);

    % Get indices of disabled parameters
    disabled_pidx = BMA.Pind(is_disabled);    
else
    disabled_pidx = [];
    BMA = [];
end

% Estimate PEBs for each dataset using the optimal DCM structure from above
% -------------------------------------------------------------------------

% Update design matrix
M.X = ones(ns,1);

PEBs = {};
for i = 1:ny
    % Load dataset
    GCM = spm_dcm_load(GCMs{i});
        
    % Disable any parameters which were pruned from the PEB.
    rE = GCM{1}.M.pE;
    rC = GCM{1}.M.pC;
    for j = 1:length(disabled_pidx)
        rC(disabled_pidx(j),disabled_pidx(j)) = 0;
    end
    
    % Apply priors at DCM level
    GCM = spm_dcm_reduce(GCM, rE, rC);
    
    % Build PEB
    M.X = ones(size(GCM,1),1);
    PEBs{i} = spm_dcm_peb(GCM,M,field);
        
end

% Compute measures
% -------------------------------------------------------------------------
d = struct();
d.precisions = [];
d.dcm_negent = [];
d.rfx_negent = [];
d.complexity = [];

for i = 1:ny       
    % Certainty about RFX
    d.rfx_negent(i) = -0.5 * spm_logdet(2 * pi * exp(1) * PEBs{i}.Ch);
    
    % Certainty about parameters
    d.dcm_negent(i) = -0.5 * spm_logdet(2 * pi * exp(1) * PEBs{i}.Cp);

    % Information gain (KL-divergence) over parameters
    qE    = PEBs{i}.Ep;    % Posterior expectations
    pE    = PEBs{i}.M.pE;  % Prior expectations
    qC    = PEBs{i}.Cp;    % Posterior covariance
    pC    = PEBs{i}.M.pC;  % Prior covariance
    
    d.complexity(i) = spm_kl_normal(qE,qC,pE,pC);
end

% Assess model discrimination
% -------------------------------------------------------------------------

% If no models are provided, enter an automatic mode where difficult to
% distinguish parameters are identified
is_automatic = isempty(models);

% Form a model space K with one row per model and one column per parameter
% where 1 = on and 0 = off.
if isempty(models)
    % Default: each parameter turned off individually
    np = length(PEBs{1}.Pnames);
    K  = ones(np, np);
    for i = 1:np
        K(i,i) = 0;
    end
elseif iscell(models) && ~isempty(models) && iscell(models{1})
    % Nested cell array of included fields
    Pnames = char(PEBs{1}.Pnames);
    
    K  = [];
    
    for m = 1:length(models)
        k = any(ismember(Pnames,models{m}),2);
        K(end+1,:) = k';
    end
    
    % Remove any models where parameters were provided but none was found
    % in the PEB.
    to_remove = ~cellfun(@isempty,models) & ~any(K,2);
    K(to_remove,:) = [];    
elseif ismatrix(models)
    % Pre-defined matrix of parameters
    K = models;    
else
    error('Unknown model specification');
end

% Calculate evidence for reduced models with each connection switched off
% Fs(m,i) is the free energy of the PEB model with parameter m disabled
% in dataset i.
Fs = []; % models x datasets
for i = 1:ny
    qE = PEBs{i}.Ep;
    qC = PEBs{i}.Cp;
    pE = PEBs{i}.M.pE;    
    pC = PEBs{i}.M.pC;
    
    for m = 1:size(K,1) 
        rE = pE;
        rC = pC;
        off = find(K(m,:) == 0);
        rC(off,off) = gamma;
        
        Fs(m,i) = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
    end
end

if is_automatic
    % Retain nr difficult to discriminate models
    i  = find(mean(Fs,2) > -3);
else
    % Retain them all
    i = 1:size(Fs,1);    
end
nr = numel(i);

% Evidence -> probability
p  = spm_softmax(Fs(i,:));

% KL divergence (information gain)
d.model_KL = sum(p.*log(p)) + log(nr);

d.model_F = Fs;
d.model_P = p;

% Get PEB parameter precisions for each dataset
% -------------------------------------------------------------------------
for i = 1:ny
    d.precisions(:,i) = diag(spm_inv(PEBs{i}.Cp));
end

if noplot
    return;
end

% Plot PEB parameters for each dataset
% -------------------------------------------------------------------------
spm_figure('GetWin','Datasets');
rows = max(ceil(ny / 2),3);
cols = 2;

for i = 1:ny
    subplot(rows,cols,i);
    spm_plot_ci(PEBs{i}.Ep,PEBs{i}.Cp);
    title(ynames{i},'FontSize',16);
    xlabel('PEB Parameter','FontSize',12);
    ylabel('Estimate','FontSize',12);
    ylim([-2 2]);
end

% Plot measures
% -------------------------------------------------------------------------
spm_figure('GetWin','Data Comparison');
spm_clf;
rows = 3; cols = 2;

subplot(rows,cols,1:2);
bar(d.precisions);
title('Parameter certainty','FontSize',16); 
xlabel('Parameter','FontSize',12); ylabel('Precision','FontSize',12);
legend(ynames,'Location','South','Orientation','horizontal','FontSize',12);

subplot(rows,cols,3);
y = d.dcm_negent - min(d.dcm_negent);
bar(y);
hold on; p = scatter(1:length(y),y); set(p,'visible','off'); hold off
title('Parameter certainty','FontSize',16);
xlabel('Dataset','FontSize',12); ylabel('Relative neg entropy (nats)','FontSize',12);

subplot(rows,cols,4);
y = d.rfx_negent - min(d.rfx_negent);
bar(y);
hold on; p = scatter(1:length(y),y); set(p,'visible','off'); hold off
title('RFX certainty','FontSize',16);
xlabel('Dataset','FontSize',12); ylabel('Relative neg entropy (nats)','FontSize',12);

subplot(rows,cols,5);
y = d.complexity-min(d.complexity);
bar(y);
hold on; p = scatter(1:length(y),y); set(p,'visible','off'); hold off
title('Information gain (parameters)','FontSize',14);
xlabel('Dataset','FontSize',12); ylabel('Relative neg entropy (nats)','FontSize',12);

subplot(rows,cols,6);
y = d.model_KL;
bar(y); 
hold on; p = scatter(1:length(y),y); set(p,'visible','off'); hold off
title(sprintf('Information gain (%d models)',nr),'FontSize',14);
xlabel('Dataset','FontSize',12); ylabel('Entropy (nats)','FontSize',12);
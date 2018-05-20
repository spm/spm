function [d,BMA,PEBs] = spm_dcm_bdc(GCMs,field,M,ynames)
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
% ynames - (Optional) {1 x ny} cell array of names for each dataset
%
% ny = number of datasets, nf = number of DCM fields
% 
% Returns:
% 
% d.precisions - [np x ny] Precision of each DCM parameter in each dataset
% d.dcm_negent - [1 x ny]  Negative entropy (certainty) of DCM parameters
% d.rfx_negent - [1 x ny]  Negative entropy (certainty) of the estimated
%                          between-subject variability
% d.complexity - [1 x ny]  Number of effective parameters in the model
% d.model_KL   - [1 x ny]  Ability to disriminate between similar models
%
% BMA        - Bayesian model average across all datasets
% PEBs       - [1 x ny] PEB for each dataset
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
 
% Peter Zeidman & Samira Kazan
% $Id: spm_dcm_bdc.m 7315 2018-05-20 10:42:51Z peter $
 
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

% Identify optimum model structure across all datasets
% -------------------------------------------------------------------------

% Perform PEB analysis
if isempty(M)
    M       = struct();
    M.alpha = 1;
    M.beta  = 16;
    M.hE    = 0;
    M.hC    = 1/16;
    M.Q     = 'single';
    M.X     = ones(nm,1);
end
PEB = spm_dcm_peb(GCM,M,field);

% Prune full model
BMA = spm_dcm_peb_bmc(PEB);

% Estimate PEBs for each dataset using the optimal DCM structure from above
% -------------------------------------------------------------------------

% Identify enabled PEB parameters
is_disabled = find(BMA.Pp < 1/1000);

% Get indices of disabled parameters
disabled_pidx = BMA.Pind(is_disabled);

% Update design matrix
M.X = ones(ns,1);

PEBs = {};
for i = 1:ny
    % Load dataset
    GCM = load(GCMs{i});
    GCM = GCM.GCM;
        
    % Build prior covariance
    rE = GCM{1}.M.pE;
    rC = GCM{1}.M.pC;
    for j = 1:length(disabled_pidx)
        rC(disabled_pidx(j),disabled_pidx(j)) = 0;
    end
    
    % Apply priors at DCM level
    GCM = spm_dcm_reduce(GCM, rE, rC);
    
    % Build PEB
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
    
    k     = rank(full(pC));
    ipC   = pinv(pC);
    
    d.complexity(i) = 0.5 * trace(ipC*qC) + (pE - qE)'*ipC*(pE - qE) - k + log(det(pC) / det(qC));
end

% Assess model discrimination
% -------------------------------------------------------------------------

% Calculate evidence for reduced models with each connection switched off
Fs = []; % models x datasets
for i = 1:ny
    qE = PEBs{i}.Ep;
    qC = PEBs{i}.Cp;
    pE = PEBs{i}.M.pE;    
    pC = PEBs{i}.M.pC;
    
    for m = 1:length(pE)    
        rE = pE;
        rC = pC;
        rC(m,m) = 0;
        
        Fs(m,i) = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
    end
end

% Retain nr difficult to discriminate models
i  = find(mean(Fs,2) > -3);
nr = numel(i);

% Evidence -> probability
p  = spm_softmax(Fs(i,:));

% KL divergence (information gain)
d.model_KL = sum(p.*log(p)) + log(nr);

% Get PEB parameter precisions for each dataset
% -------------------------------------------------------------------------
for i = 1:ny
    d.precisions(:,i) = diag(spm_inv(PEBs{i}.Cp));
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
bar(d.dcm_negent - min(d.dcm_negent));
title('Parameter certainty','FontSize',16);
xlabel('Dataset','FontSize',12); ylabel('Relative neg entropy (nats)','FontSize',12);

subplot(rows,cols,4);
bar(d.rfx_negent - min(d.rfx_negent));
title('RFX certainty','FontSize',16);
xlabel('Dataset','FontSize',12); ylabel('Relative neg entropy (nats)','FontSize',12);

subplot(rows,cols,5);
bar(d.complexity-min(d.complexity));
title('Information gain (parameters)','FontSize',14);
xlabel('Dataset','FontSize',12); ylabel('Relative neg entropy (nats)','FontSize',12);

subplot(rows,cols,6);
bar(d.model_KL-min(d.model_KL));
title(sprintf('Information gain (%d models)',nr),'FontSize',14);
xlabel('Dataset','FontSize',12);
xlabel('Dataset','FontSize',12); ylabel('Relative entropy (nats)','FontSize',12);
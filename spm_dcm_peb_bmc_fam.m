function [BMA,fam] = spm_dcm_peb_bmc_fam(BMA,BMR,families,bma_option)
% Bayesian model selection and averaging over families of PEB models
%
% BMA - Bayesian model average (see spm_dcm_peb_bmc)
% 
%            BMA.K       - model space [models x parameters]
%            BMA.Pnames  - parameter names
%            BMA.Ep      - averaged parameters (for inferring #covariates)
%
% BMR - model space (see spm_dcm_peb_bmc)
%
%            BMR{i,j}    - model i of commonalities and j of differences
%            BMR{i,j}.Ep - expectations of second level parameters
%            BMR{i,j}.Cp - covariance of second level parameters
%            BMR{i,j}.F  - free energy relative to full model
%
% families - [1 x Nm] vector of family membership where families(i)=x,
%            for model i and family x. For example, families=[1 1 2]
%            means that models 1 and 2 are in family 1 and model 3 is in
%            family 2.
%
% bma_option - String specifying option for Bayesian Model Averaging
%
%            'ALL'     - average over all models (under family priors)
%            'WINNING' - average models in the best family only
%            'NONE'    - don't average
%
% Returns:
%
% BMA - Bayesian Model Average over models (specified by bma_option)
%
% fam - Bayesian model selection results:
%   
%           .model.post   - Posterior probability of each model
%           .model.prior  - Prior probability of each model
%           .family.post  - Posterior probability of each family
%           .family.prior - Prior probability of each family
%           .famdef       - Input vector of family assignments
%__________________________________________________________________________
% Copyright (C) 2010-2016 Wellcome Trust Centre for Neuroimaging

% Peter Zeidman
% $Id: spm_dcm_peb_bmc_fam.m 6946 2016-11-23 15:26:29Z peter $

% Model space [models x parameters]
K      = BMA.K;
Nm     = size(K,1);

% Number of models for between-subjects effects
Nmx    = size(BMR,2);

% Validate
if nargin < 4
    bma_option = 'ALL';
end

if ~iscell(BMR)
    error('BMR should be a cell array of parameters and evidences');
end

if ~isvector(families) || length(families) ~= Nm
    error('Families should be a vector with one element per model');
end

% Unpack
%--------------------------------------------------------------------------
for i = 1:Nm
    for j = 1:Nmx
        F(i,j) = BMR{i,j}.F;        
    end    
end

Pnames = BMA.Pnames;                           % Parameter names
M      = BMA.M;                                % Original PEB model
Nx     = length(BMA.Ep) / length(BMA.Pnames);  % Number of covariates
nK     = length(unique(families));             % Number of families   

% Set model priors given family priors
%--------------------------------------------------------------------------

% Likelihood
L = exp(F(:) - max(F(:)));

% Size of each family
sK = zeros(1,nK);
for i = 1:nK
    sK(i) = sum(families==i);
end

% Family prior
pF = repmat(1/nK, 1, nK);

% Model prior (commonalities)
pMw = 1 ./ (sK(families) * nK);

% Model prior (differences)
if Nmx > 1
    pMb = pMw;
else
    pMb = 1;
end

% Joint model prior over commonalities and differences
pMj = pMw' * pMb;

% Calculate posteriors
%--------------------------------------------------------------------------

% Model posterior
pMj = pMj(:);
Pm = L .* pMj;
Pm = Pm ./ sum(Pm);
Pm = reshape(Pm,Nm,Nmx);

% Group differences families (fx)
if Nmx > 1
    fx = families;
    k  = nK;
else
    fx = 1;
    k  = 1;
end

% Family posterior
Pf    = zeros(nK,k);
for i = 1:nK   
    for j = 1:k
        p = Pm(families == i, fx == j);
        Pf(i,j) = sum(p(:));
    end
end

% Compute BMA
%--------------------------------------------------------------------------
bma_option = upper(bma_option);

occ = 8; % Occam's window

switch bma_option
    case 'ALL'
        i   = F(:) > max(F(:) - occ);
        m   = reshape(i,Nm,Nmx);
        BMA = spm_dcm_bma(BMR(i)');
    case 'WINNING'
        % Identify winning family
        [wF,i]    = max(Pf(:));
        [wFi,wFj] = ind2sub([nK nK],i);
        
        % Identify models
        i = (families == wFi);
        j = (fx == wFj);
                
        % Matrix of included models for BMA
        m = zeros(Nm,Nmx);
        m(i,j) = 1;
        m      = m .* (F > max(F(:) - occ));
        
        % Compute BMA
        BMA = spm_dcm_bma(BMR(m(:)==1)');        
    case 'NONE'
        BMA = [];
    otherwise
        error('Unknown BMA option');
end

% Inference over models (commonalities and differences)
%--------------------------------------------------------------------------
P     = Pm;
P1    = sum(P,2);
P2    = sum(P,1);

% Inference over parameters (present vs absent)
%--------------------------------------------------------------------------
k     = find(any(~K));
Nb    = length(k);
Kname = Pnames(k);
for i = 1:Nb
    Pw(1,i) = mean(sum(P( ~K(:,k(i)),:),2));
    Pw(2,i) = mean(sum(P(~~K(:,k(i)),:),2));
    
    if Nx > 1
        Px(1,i) = mean(sum(P(:, ~K(:,k(i))),1));
        Px(2,i) = mean(sum(P(:,~~K(:,k(i))),1));
    else
        Px(1,i) = 1;
        Px(2,i) = 0;
    end
    
end
Pw    = Pw(2,:)./sum(Pw,1);
Px    = Px(2,:)./sum(Px,1);

% Pack
%--------------------------------------------------------------------------
BMA.F      = F;
BMA.P      = P;
BMA.Px     = Px;
BMA.Pw     = Pw;
BMA.M      = M;
BMA.K      = K;
BMA.Pnames = Pnames;
BMA.Kname  = Kname;

fam = struct();
fam.model.post   = Pm;
fam.model.prior  = reshape(pMj,Nm,Nmx);
fam.model.Pw     = P1';
fam.model.Px     = P2;
fam.family.post  = Pf;
fam.family.prior = pF;
fam.famdef       = families;

% Plot
%--------------------------------------------------------------------------
if spm_get_defaults('cmdline'), return; end

spm_figure('GetWin','Family analysis');
clf;

% Family posterior
subplot(3,2,1)
if Nmx > 1
    imagesc(Pf)
    xlabel('Family (differences)','FontSize',12)
    ylabel('Family (commonalities)','FontSize',12)    
    set(gca,'XTick',1:nK,'YTick',1:nK);
else
    bar(Pf);
    xlabel('Family','FontSize',12);
    ylabel('Probability');
end
title('Family posterior','FontSize',16)
axis square

% Model posterior
subplot(3,2,2)
if Nmx > 1
    imagesc(Pm)
    xlabel('Model (differences)','FontSize',12)
    ylabel('Model (commonalities)','FontSize',12)
else
    bar(Pm);
    xlabel('Model','FontSize',12);
end
title('Model posterior | families','FontSize',16)
axis square

sp = 3; % Subplot counter

if Nmx > 1
    % Commonalities
    subplot(3,2,sp)
    [mm,i] = max(P1); bar(P1),
    text(i - 1/4,mm/2,sprintf('%-2.0f%%',mm*100),'Color','w','FontSize',8)
    title('Commonalities','FontSize',16)
    xlabel('Model','FontSize',12)
    ylabel('Probability','FontSize',12)
    axis([0 (Nm + 1) 0 1]), axis square
    sp = sp + 1;

    % Differences
    subplot(3,2,sp)
    [mm,i] = max(P2); bar(P2),
    text(i - 1/4,mm/2,sprintf('%-2.0f%%',mm*100),'Color','w','FontSize',8)
    title('Differences','FontSize',16)
    xlabel('Model','FontSize',12)
    ylabel('Probability','FontSize',12)
    axis([0 (Nm + 1) 0 1]), axis square
    sp = sp + 1;
end

% Family allocations
subplot(3,2,sp), imagesc(families)
title('Families','FontSize',16)
xlabel('Model','FontSize',12)
set(gca,'YTick',[]);
axis square
sp = sp + 1;

% BMA inclusion
if ~strcmp(bma_option,'NONE')
    subplot(3,2,sp), imagesc(m)
    title('Models in BMA','FontSize',16)
    xlabel('Model (differences)','FontSize',12)
    ylabel('Model (commonalities)','FontSize',12)
    
    if Nmx == 1
        set(gca,'XTick',1);
    end
end
axis square
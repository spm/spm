function [BMA,BMR] = spm_dcm_peb_bmc(PEB,models)
% hierarchical (PEB) model comparison and averaging (2nd level)
% FORMAT [BMA] = spm_dcm_peb_bmc(PEB,models)
% FORMAT [BMA] = spm_dcm_peb_bmc(PEB)
%
% PEB -  between subject (second level) effects (from spm_dcm_peb)
% ------------------------------------------------------------
%     PEB.Snames - string array of Ns first level model names
%     PEB.Pnames - string array of Np parameters of interest
%
%     PEB.M.X  -   second level (between subject) design matrix
%     PEB.M.W  -   second level (within  subject) design matrix
%     PEB.M.Q  -   precision components of second level random effects
%     PEB.M.pE -   prior expectation of second level parameters
%     PEB.M.pC -   prior covariance  of second level parameters
%     PEB.Ep   -   posterior expectation of second level parameters
%     PEB.Cp   -   posterior covariance  of second level parameters
%
% models - field in DCM.Ep to compare For the first two group effects
%          or logical (Nm x Np) matrix of Nm (parametric) model space
%          or an array of DCMs specifying Nm (parametric) model space
%
% if models are not specified, all combinations of second level parameters
% will be tested.
%
% BMA    - DCM structure of Bayesian model average
% -------------------------------------------------------------------------
%     BMA.Snames - string array of first level model names
%     BMA.Pnames - string array of parameters of interest
%
%     BMA.Ep   - BMA expectation of second level parameters
%     BMA.Cp   - BMA   variances of second level parameters
%     BMA.M    - second level model
%
%     BMA.F    - free energy over model space
%     BMA.P    - posterior probability over models
%     BMA.Px   - posterior probability over parameters (differences)
%     BMA.Pw   - posterior probability over parameters (common)
%     BMA.K    - model space
%
% or for automatic model search, see spm_dcm_bmr_all.m (output: DCM)
%
% BMR   - Parameters and evidence of reduced models which produced the BMA
% -------------------------------------------------------------------------
%     BMR{i,j}    - model i of commonalities and j of group differences
%     BMR{i,j}.Ep - expectations of second level parameters
%     BMR{i,j}.Cp - covariance of second level parameters
%     BMR{i,j}.F  - free energy relative to full model
%
% or for automatic model search:
%
%     BMR.name   - parameter names
%     BMR.F      - free energy relative to full model
%     BMR.P      - and posterior (model) probabilities
%     BMR.K      - [models x parameters] model space (1 = off, 2 = on)
%     BMR.bma{i} - Model i which contributed to the BMA (Ep,Cp,F)
%
%--------------------------------------------------------------------------
% This routine performs Bayesian model comparison and averaging of second
% level or hierarchical (PEB) models. The model space is defined either
% in terms of fields (e.g. 'A' or 'B') or as a logical matrix, with one row
% per model and a column per parameter (in PEB.Pnames). This induces
% a joint model space over parameters and group effects at the second level
% (encoded by the design matrix, X). Using Bayesian model reduction, this
% joint model space is scored over the specified models at the first level
% (for the constant terms modelling effects that are common to all
% subjects) and combinations of group effects (modelling between
% subject differences).
%
% If there is only a group effect (and no between subject differences) this
% reduces to a search over different models of the group mean.
%
% Given the model space one can then compute the posterior probability
% of various combinations of group effects over different parameters. Of
% particular interest are (i) the posterior probabilities over the
% the first two group effects in the design matrix and the posterior
% probability of models with and without each parameter, for the common
% (first) and subject-specific (second) group affects (returned in BMA.P,
% BMA.Pw and BMA.Px respectively. The Bayesian model averages of the second
% level parameters and can be found in BMA.Ep and BMA.Cp.
%
% If models are not specified, all combinations of individual
% parameters over all group effects will be considered and the ensuing
% Bayesian model reduction reported for each effect in the design matrix.
%
% NB for EEG models the absence of a connection means it is equal to its
% prior mesn, not that is is zero.
%
% see also: spm_dcm_peb.m and spm_dcm_bmr
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_peb_bmc.m 7081 2017-05-27 19:36:09Z karl $

% checks
%--------------------------------------------------------------------------
if nargin < 1 || isempty(PEB) || length(PEB) > 1
    error('Please provide a single PEB model');
end

% (greedy) search over all combinations of second level parameters
%==========================================================================
if nargin < 2
    
    % greedy search and (second level) Bayesian model average
    %----------------------------------------------------------------------
    [BMA,BMR,bma]  = spm_dcm_bmr_all(PEB);
       
    % plot posteriors over parameters
    %======================================================================
    if spm_get_defaults('cmdline'), return; end
    
    spm_figure('Getwin','BMC'); clf
    
    Np    = numel(PEB.Pind);
    Nx    = min(3,size(PEB.M.X,2));
    
    if isfield(PEB,'Xnames')
        str   = PEB.Xnames;
    else
        str   = {'Group mean','1st group effect','2nd group effect'};
    end
        
    for i = 1:Nx
        
        j = (1:Np)' + (i - 1)*Np;
        
        % posterior density over parameters
        %------------------------------------------------------------------
        subplot(3,Nx,0 + i), spm_plot_ci(PEB.Ep(j),PEB.Cp(j,j))
        title(str{i},'FontSize',16)
        xlabel('Parameter','FontSize',12)
        ylabel('Effect size','FontSize',12)
        axis square, a = axis;
        if var(PEB.Ep(j)) < 1e-6
            a(3:4) = [-1 1];
            axis(a)
        end
        
        
        % posterior density over parameters
        %------------------------------------------------------------------
        subplot(3,Nx,Nx + i), spm_plot_ci(bma.Ep(j),bma.Cp(j))
        title('Reduced','FontSize',16)
        xlabel('Parameter','FontSize',12)
        ylabel('Effect size','FontSize',12)
        axis square, axis(a);

        % posterior density over parameters
        %------------------------------------------------------------------
        subplot(3,Nx,Nx + Nx + i)
        if Np > 1
            bar(diag(BMA.Pp(j)),Np)
        else
            bar(BMA.Pp(j))
        end
        title('Posterior','FontSize',16)
        xlabel('Parameter','FontSize',12)
        ylabel('Probability','FontSize',12)
        axis([0 (Np + 1) 0 1]), axis square
    
    end
    legend(BMA.Pnames)
    return
end

% otherwise search a smaller joint space of first and second group effects
%==========================================================================


% number of parameters and effects
%--------------------------------------------------------------------------
[Np,Nx]   = size(PEB.Ep);

if ischar(models)
    
    % compare all combinations of field in 'models'
    %----------------------------------------------------------------------
    Pnames = char(PEB.Pnames);
    k      = any(ismember(Pnames,models),2);
    K      = ones(2^sum(k),Np);
    K(:,k) = spm_perm_mtx(sum(k));
    
elseif iscell(models)
    
    % (RFX) BMA – define the model space in terms of a matrix
    %----------------------------------------------------------------------
    Nm    = length(models);
    Np    = length(PEB.Pind);
    K     = ones(Nm,Np);
    for i = 1:Nm
        k      = spm_find_pC(models{i});
        j      = find(~ismember(PEB.Pind,k));
        K(i,j) = 0;
    end
    
else
    
    % model space in defined in terms of a matrix
    %----------------------------------------------------------------------
    K     = models;
    
end
[k,i]     = unique(K,'rows');
K         = K(sort(i),:);
[Nm,Np]   = size(K);


% check number of models
%--------------------------------------------------------------------------
i = find(any(~K),1);
if isempty(i)
    warndlg('your model space is empty')
    return
end


%-score models with log-evidences
%==========================================================================
fprintf('BMC:     ')

% Get priors and posteriors - of first and second order parameters
%--------------------------------------------------------------------------
qE    = spm_vec(PEB.Ep);
qC    = PEB.Cp;
pE    = spm_vec(PEB.M.pE);
pC    = PEB.M.pC;
for i = 1:Nm
    
    if Nx > 1
        
        % model comparison over common (constant) and group effects
        %------------------------------------------------------------------
        for j = 1:Nm
            
            % reduced prior
            %--------------------------------------------------------------
            k   = [K(i,:) K(j,:) ones(1,(Nx - 2)*Np)];
            R   = diag(k);
            rE  = R*pE;
            rC  = R*pC*R;
            
            % Bayesian model reduction (of second level)
            %--------------------------------------------------------------
            [F, sE, sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
            BMR{i,j}.Ep = sE;
            BMR{i,j}.Cp = sC;
            BMR{i,j}.F  = F;
            G(i,j)      = F;
            
            % report progress
            %--------------------------------------------------------------
            fprintf('\b\b\b\b%-3.0f%%',100*((i - 1)*Nm + j)/(Nm*Nm))
            
        end
        
    else
        
        % otherwise, reduced prior over group mean
        %------------------------------------------------------------------
        k   = K(i,:);
        R   = diag(k);
        
        rE  = R*pE;
        rC  = R*pC*R;
        
        % Bayesian model reduction (of second level)
        %------------------------------------------------------------------
        [F, sE, sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);
        BMR{i,1}.Ep = sE;
        BMR{i,1}.Cp = sC;
        BMR{i,1}.F  = F;
        G(i,1)      = F;
    end
end

% family wise inference over models and parameters
%==========================================================================
P     = G;
P(:)  = exp(P(:) - max(P(:)));
P(:)  = P/sum(P(:));

% family wise inference over parameters (present an absent)
%--------------------------------------------------------------------------
k     = find(any(~K));
Nb    = length(k);
Kname = PEB.Pnames(k);
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

% family wise inference over mmodels (commonalities and differences)
%--------------------------------------------------------------------------
P1    = sum(P,2);
P2    = sum(P,1);

%-hierarchical inversion using optimised second level priors
%==========================================================================

% Bayesian model averaging (with an Occam's window of eight)
%--------------------------------------------------------------------------
i     = G(:) > max(G(:) - 8);
BMA   = spm_dcm_bma(BMR(i)');

% assemble BMA output structure
%--------------------------------------------------------------------------
BMA.Snames = PEB.Snames;
BMA.Pnames = PEB.Pnames;
BMA.Pind   = PEB.Pind;
BMA.Kname  = Kname;
try BMA.Xnames = PEB.Xnames; catch, BMA.Xnames = {}; end

BMA.F     = G;
BMA.P     = P;
BMA.Px    = Px;
BMA.Pw    = Pw;
BMA.M     = PEB.M;
BMA.K     = K;

% add posterior precisions from PEB (not averaged)
%--------------------------------------------------------------------------
BMA.Ce    = PEB.Ce;
BMA.Ch    = PEB.Ch;
BMA.Eh    = PEB.Eh;

% Show results
%==========================================================================

% get name of covariate 2 (differences)
%--------------------------------------------------------------------------
if Nx > 1 && isfield(PEB,'Xnames') && ~isempty(PEB.Xnames)
    xname = PEB.Xnames{2};
else
    xname = 'Differences';
end

if spm_get_defaults('cmdline'), return; end

spm_figure('Getwin','BMC'); clf

subplot(3,2,1), imagesc(K')
title('Model space','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Parameter','FontSize',12)
set(gca,'YTick',1:Np,'YTickLabel',PEB.Pnames)
axis square

subplot(3,2,3)
[m,i] = max(P1); bar(P1),
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Commonalities','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (Nm + 1) 0 1]), axis square

subplot(3,2,5), bar(diag(Pw),length(Pw));
title('Commonalities','FontSize',16)
xlabel('Parameter','FontSize',12)
ylabel('Parameter probability','FontSize',12)
axis([0 (Nb + 1) 0 1]), axis square
legend(Kname)

subplot(3,2,2), imagesc(P)
title('Posterior probabilities','FontSize',16)
xlabel('Model (differences)','FontSize',12)
ylabel('Model (commonalities)','FontSize',12)
axis square

if Nx < 2
    
    % posterior density over parameters
    %----------------------------------------------------------------------
    subplot(3,2,4), spm_plot_ci(PEB.Ep,PEB.Cp)
    title('MAP (full)','FontSize',16)
    xlabel('Parameter','FontSize',12)
    ylabel('Effect size','FontSize',12)
    axis square
    
    % posterior density over parameters
    %----------------------------------------------------------------------
    subplot(3,2,6), spm_plot_ci(BMA.Ep,BMA.Cp)
    title('BMA (reduced)','FontSize',16)
    xlabel('Parameter','FontSize',12)
    ylabel('Effect size','FontSize',12)
    axis square
    
else
    
    % inference over group effects
    %----------------------------------------------------------------------
    subplot(3,2,4)
    [m,i] = max(P2); bar(P2),
    text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
    title(xname,'FontSize',16)
    xlabel('Model','FontSize',12)
    ylabel('Probability','FontSize',12)
    axis([0 (Nm + 1) 0 1]), axis square
    
    subplot(3,2,6), bar(diag(Px),length(Px))
    title(xname,'FontSize',16)
    xlabel('Parameter','FontSize',12)
    ylabel('Parameter probability','FontSize',12)
    axis([0 (Nb + 1) 0 1]), axis square
    
end




